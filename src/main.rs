use concrete_core::commons::crypto::encoding::PlaintextList;

use panacea::context::*;
use panacea::num_types::*;
use panacea::params::*;
use panacea::rlwe::*;

fn main() {
    let mut ctx = Context::new(TFHEParameters::default());
    let ptxt_expected = ctx.gen_binary_pt();

    let sk = ctx.gen_rlwe_sk();
    let mut ct = RLWECiphertext::allocate(ctx.poly_size);
    sk.encode_encrypt_rlwe(&mut ct, &ptxt_expected, &mut ctx);

    let mut ptxt_actual = PlaintextList::allocate(Scalar::zero(), ctx.plaintext_count());
    sk.decrypt_decode_rlwe(&mut ptxt_actual, &ct, &ctx);

    println!("ok? {}", ptxt_actual == ptxt_expected);
}
