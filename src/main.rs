use concrete_core::commons::crypto::encoding::PlaintextList;
use concrete_core::commons::crypto::glwe::GlweCiphertext;
use concrete_core::commons::math::polynomial::{MonomialDegree};
use concrete_core::commons::math::tensor::{AsMutTensor, AsRefTensor};
use panacea::context::*;
use panacea::num_types::*;
use panacea::params::*;
use panacea::rlwe::*;
use panacea::rgsw::*;
use std::collections::HashMap;
use std::time::Instant;
use std::io;

struct Queue<T> {
    queue: Vec<T>,
}

impl<T> Queue<T> {
    fn new() -> Self {
        Queue { queue: Vec::new() }
    }

    fn enqueue(&mut self, item: T) {
        self.queue.push(item)
    }

    fn dequeue(&mut self) -> T {
        self.queue.remove(0)
    }
}
//------------------------------------------------------------------------------------------------------------------------------------
// Function for trivial encryption

fn encode_trivial_encrypt_rlwe(pt: &PlaintextList<Vec<Scalar>>, ctx: &mut Context)-> RLWECiphertext{
    let mut binary_encoded = pt.clone();
        ctx.codec
            .poly_encode(&mut binary_encoded.as_mut_polynomial());

    let mut ct = RLWECiphertext::allocate(ctx.poly_size);
    ct.0.fill_with_trivial_encryption(&binary_encoded);    
    return ct;
}

//------------------------------------------------------------------------------------------------------------------------------------
//Function for generating the client's query. 
// If query index = SUM\i=0^n-1(e_i*2^i), then the function outputs RLWE(SUM\i=0^n-1(q/N*B_g^j)*e_i*x^i)) for all j=1(1)l    -  N = polynomial degree
fn query_gen(x: u64, n: usize, ctx: &mut Context, sk: &RLWESecretKey) -> Vec<RLWECiphertext>{			//n = log_2(database size)
    
    let base: i32 = 2;
    let mut dummy_x : u64 = x;   
    let l = ctx.level_count.0;
    let b_g = ctx.base_log.0;
    let p = ctx.poly_size.0.ilog2() as usize;
   // println!("polynomial aize log - {}",p);
    let mut query: Vec<RLWECiphertext> = Vec::with_capacity(l);
    
    for i in 0..l {
        query.insert(i, RLWECiphertext::allocate(ctx.poly_size));
    } 
    
     // For printing the noise - for cross-checking
    //println!("Noise for each query polynomial");
    //println!("---------------------------------");
    
    for i in 0..l{
        dummy_x = x;
        let mut temp = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
        let shift: usize = (Scalar::BITS as usize) - p - (b_g * (i+1));
        
        for j in 0..n{
            *temp
            .as_mut_polynomial()
            .get_mut_monomial(MonomialDegree(j))
            .get_mut_coefficient() = (dummy_x & 1) * (1<<shift);
            dummy_x = dummy_x>>1;
        }
        //println!("{}>>{:?}",i, temp);
      
        sk.encrypt_rlwe(&mut query[i], &temp, ctx.std, &mut ctx.encryption_generator);
       
       /* Printing the decryption of the encrypted query for comparison purpose 
        let mut decrypted = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
        sk.decrypt_rlwe(&mut decrypted, &query[i]);
        println!("{}>>{:?}",i,query[i]);
        println!("{}>>{:?}",i, decrypted);
        println!("For {}th query:: Result Matching? {}	||  {:?}",i, decrypted==temp, compute_noise(&sk, &query[i], &temp));
        println!("---------------------------------");
        
        //println!("For {}th query: {:?}",i, compute_noise(&sk, &query[i], &temp));
      */
    }
    return query;    
}

//--------------------------------------------------------------------------------------------------------------------------------
//From the l RLWEciphertexts, the following function extracts RLWE((q/B_g^j)*e_i), for all i=0(1)n-1 , for all j=1(1)l by rotating the rlwe ciphertexts and using trace function on them. 
// trace1(RLWE(\sum_i` `a_i` X^i)) = `RLWE(N*a_0`)   - N polynomial degree
fn conversion_trace(query: &mut Vec<RLWECiphertext>, n: usize, ctx: &mut Context, ksk_map: &HashMap<usize, FourierRLWEKeyswitchKey>)-> Vec<RLWECiphertext>{

    let l = ctx.level_count.0;
    let mut out: Vec<RLWECiphertext> = Vec::with_capacity(l*n);
    for i in 0..(l*n) {
        out.insert(i, RLWECiphertext::allocate(ctx.poly_size));
    } 
    let mut temp_query: Vec<RLWECiphertext> = Vec::with_capacity(l);
    let mut temp_query = query.clone();
    
    for i in 0..n{
        for j in 0..l{
            temp_query = query.clone();
            temp_query[j].update_with_monomial_div(MonomialDegree(i));
            //out[i*l+j] = temp_query[j].trace1(&ksk_map);  
            temp_query[j].trace1_fourier(&mut out[i*l+j], &ksk_map);  
        }
    }     
    return out; 
}

//---------------------------------------------------------------------------------------------------------------------------------
//The above function gives l*n many RLWE ciphertexts from l rlwe obtained from client. The following function takes l RLWE ciphertexts corresponding to any bit e_i and returns RGSW(e_i).

fn conversion_expand1(cs_prime: &Vec<RLWECiphertext>, neg_s: &RGSWCiphertext, ctx: &Context) -> RGSWCiphertext {
    let mut buf = FftBuffer::new(cs_prime[0].polynomial_size());
    
    let mut out = RGSWCiphertext::allocate(ctx.poly_size, ctx.base_log, ctx.level_count);
    for (i, mut c) in out.0.as_mut_glwe_list().ciphertext_iter_mut().enumerate() {
        let k = i / 2;
        if i % 2 == 0 {
            neg_s.external_product_with_buf_glwe(&mut c, &cs_prime[k], &mut buf);
        } else {
            c.as_mut_tensor().fill_with_copy(cs_prime[k].0.as_tensor());
        }
    }
    return out;
}

//---------------------------------------------------------------------------------------------------------------------------------
//Using the avobe function, it outputs n RGSWs one corresponding to each bit from the l RLWEs obtained from the client
fn conversion_expand2(c:&Vec<RLWECiphertext>, n: usize, neg_s: &RGSWCiphertext, ctx: &Context)-> Vec<RGSWCiphertext>{
    let l = ctx.level_count.0;
    let mut id_ctx: Vec<RGSWCiphertext> = Vec::with_capacity(n);
    for i in 0..n{
        id_ctx.insert(i as usize, RGSWCiphertext::allocate(ctx.poly_size, ctx.base_log, ctx.level_count));
    }
    
    for i in 0..n{
        id_ctx[i] = conversion_expand1(&c[i*l..(i+1)*l].to_vec(), &neg_s, &ctx);
    } 
    return id_ctx;
}

//---------------------------------------------------------------------------------------------------------------------------------
// SHECS-PIR Evaluation Function

fn evaluation(id_ctx: &Vec<RGSWCiphertext>, enc_db: &Vec<RLWECiphertext>, N: usize, sk: &RLWESecretKey, ctx: &Context) -> RLWECiphertext /*PlaintextList<Vec<u64>>*/{
 
    //Preparing the queue - 1st level of evaluation
    let mut db_queue: Queue<RLWECiphertext> = Queue::new();
    let mut temp = RLWECiphertext::allocate(ctx.poly_size);   
    let n: usize = N.ilog2() as usize;
    let base: i32 = 2;
           
    for i in 0..(N/2){
        id_ctx[0].cmux(&mut temp, &enc_db[2*i], &enc_db[2*i +1]);
        db_queue.enqueue(temp.clone());
    }
   
    for i in 1..n{
        for _j in 0..(N/base.pow((i+1) as u32) as usize){
            let ct_0 = db_queue.dequeue();
            let ct_1 = db_queue.dequeue();
            id_ctx[i as usize].cmux(&mut temp, &ct_0, &ct_1);
            db_queue.enqueue(temp.clone());
        } 
    }
        
    let result_enc = db_queue.dequeue();
    return result_enc;
    /*let mut result = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    sk.decrypt_decode_rlwe(&mut result, &result_enc, &ctx);
    println!("Noise: {}", compute_noise(&sk, &result_enc, &result));*/
    //return result;
}

//----------------------------------------------------------------------------------------------------------------------------------
// Function to check the generated RGSWs
fn check(gsw_ct: &RGSWCiphertext, sk: &RLWESecretKey, ctx: &mut Context){
    
    let mut lwe_ct = RLWECiphertext::allocate(ctx.poly_size);
    let temp = ctx.gen_unit_pt();
    let mut test_ct = RLWECiphertext::allocate(ctx.poly_size);
    sk.encode_encrypt_rlwe(&mut test_ct, &temp, ctx);
    
    /*let mut pt1 = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    sk.decrypt_rlwe(&mut pt1,&gsw_ct.get_nth_row(1));
    println!("1st row: {:?}", pt1);*/   							//printing 1st row
    
    gsw_ct.external_product(&mut lwe_ct, &mut test_ct);
    let mut pt = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    sk.decrypt_decode_rlwe(&mut pt, &lwe_ct, ctx);
    println!("{:?}",pt);
}

//----------------------------------------------------------------------------------------------------------------------------------
fn main() {
    
    let base: i32 = 2;
    
    println!("Enter the value of log(DatabaseSize): ");
    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Not a valid string");
    let input_num: u32 = input.trim().parse().expect("Not a valid string");
    
    println!("Enter the value of the query index: ");
    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Not a valid string");
    let id: u64 = input.trim().parse().expect("Not a valid string");		// query index - input from  client side
    
    let N: usize= base.pow(input_num).try_into().unwrap();			//N = Database size ~ 1024
    let n: usize = input_num as usize;						// n = log_2 (N) ~ 10
   
    //let mut ctx = Context::new(TFHEParameters::default());
    let mut ctx = Context::new(TFHEParameters {
            standard_deviation: -55.0,
            polynomial_size: 2048,
            base_log: 7,		//B_g
            level_count: 7,		//l
            key_switch_base_log: 7,
            key_switch_level_count: 7,
            negs_base_log: 7,
            negs_level_count: 7,
            plaintext_modulus: 1 << 8,
            secure_seed: true,
        });
    let l = ctx.level_count.0;
    let b_g = ctx.base_log.0;
    
    // creating the parameters
    let sk = ctx.gen_rlwe_sk();
    let neg_sk_ct = sk.neg_gsw(&mut ctx);
    //let ksk_map = gen_all_subs_ksk(&sk, &mut ctx);
    let ksk_map = gen_all_subs_ksk_fourier(&sk, &mut ctx);
    
    // Creating Database
    println!("Creating Database");
    let mut db: Vec<PlaintextList<Vec<u64>>> = Vec::with_capacity(N);
    for i in 0..N {
        db.insert(i, ctx.gen_binary_pt());
    }
    
    println!("Encrypting database");
    let mut enc_db: Vec<RLWECiphertext> = Vec::with_capacity(N);
    for i in 0..N {
        enc_db.insert(i, RLWECiphertext::allocate(ctx.poly_size));
    } 
    for i in 0..N {
        //sk.encode_encrypt_rlwe(&mut enc_db[i], &db[i], &mut ctx);		// Encrypting database
        //enc_db[i].0.fill_with_trivial_encryption(&db[i]);			-- trivial encryption - not working !!!
        enc_db[i] = encode_trivial_encrypt_rlwe(&db[i], &mut ctx);
    } 
    
    println!("1. Database created");
       
    //let id: u64 = 23;										
    
    let start = Instant::now();
    let mut query: Vec<RLWECiphertext> = query_gen(id, n, &mut ctx, &sk);			//generating query rlwe ciphertexts - QUERY PACKING
    let duration1 = start.elapsed();
    println!("2. Query Generated. Time taken = {:?}",duration1);
    
    let start = Instant::now();
    let coeff = conversion_trace(&mut query, n, &mut ctx, &ksk_map);				// trace
    //println!("3. trace done.");
    
    /*//noise checking after trace
    let mut pt = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    let mut dummy_x = id;
    //let mut decrypted = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    
    println!("---------------------------------");
    println!("Noise after calculating trace");
    
    for i in 0..n{
        for j in 0..l{
            
            let shift: usize = (Scalar::BITS as usize) - 20 - (b_g * (j+1));
            *pt
            .as_mut_polynomial()
            .get_mut_monomial(MonomialDegree(0))
            .get_mut_coefficient() = (dummy_x & 1) * (1<<shift);
             
            // sk.decrypt_rlwe(&mut decrypted, &coeff[i*l+j]);
            // println!("For e_{}/(B_g)^{} :: Result Matching? {}	|| {:?} || Noise : {:?}",i, j, decrypted==pt ,decrypted, compute_noise(&sk, &coeff[i*l + j], &pt));         
            println!("For e_{}/(B_g)^{} : Noise : {:?}",i, j, compute_noise(&sk, &coeff[i*l + j], &pt));           
        }
        dummy_x = dummy_x >>1;
    }
    println!("---------------------------------");
    */
    
    let id_ctx = conversion_expand2(&coeff, n, &neg_sk_ct, &ctx);					//expansion to rgsw
    let duration2 = start.elapsed();
    println!("4. Conversion complete");
    println!("Time taken for query unpacking = {:?}", duration2);
    
     //Noise of RGSWs
    //println!("---------------------------------");
    //println!("Noise of each RGSW:");
    
    let mut max = 0.0;
    for i in 0..n{
        if compute_noise_rgsw1(&sk, &id_ctx[i], &ctx)>max{
            max = compute_noise_rgsw1(&sk, &id_ctx[i], &ctx);
        }
        //println!("Noise of e_{}: {:?}",i, compute_noise_rgsw1(&sk, &id_ctx[i], &ctx));
    }  
    println!("Maximum noise of RGSW coefficients: {}", max);
    
    /*
    for i in 0..n{
        println!("------------------ {}th bit -----------------",i);
    	check(&id_ctx[i], &sk, &mut ctx);
    }*/
    
  
    
    let start = Instant::now();
    let result_enc = evaluation(&id_ctx, &enc_db, N, &sk, &ctx);			//Query Evaluation
    let duration3 = start.elapsed();
    let mut result = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    sk.decrypt_decode_rlwe(&mut result, &result_enc, &ctx);
    		    
    println!("5. Evaluation done. Time taken for evaluation = {:?}", duration3);
    println!("Noise: {}", compute_noise(&sk, &result_enc, &result));				
    println!("ok? {}", result == db[id as usize]);   
        
}
   
   
   
   
   
   
   
   
   
   
