use crate::{ec_isogeny, f2elm, f2elm_t, fpx, PointProj, F2ELM, NWORDS_FIELD, U64};
use std::cell::RefCell;
use std::rc::Rc;
use zeroize::Zeroize;

pub fn init_basis(gen: &[u64], XP: &mut f2elm_t, XQ: &mut f2elm_t, XR: &mut f2elm_t) {
    XP[0].clone_from_slice(&gen[..NWORDS_FIELD]);
    XP[1].clone_from_slice(&gen[NWORDS_FIELD..NWORDS_FIELD * 2]);
    XQ[0].clone_from_slice(&gen[NWORDS_FIELD * 2..NWORDS_FIELD * 3]);
    XQ[1].clone_from_slice(&gen[NWORDS_FIELD * 3..NWORDS_FIELD * 4]);
    XR[0].clone_from_slice(&gen[NWORDS_FIELD * 4..NWORDS_FIELD * 5]);
    XR[1].clone_from_slice(&gen[NWORDS_FIELD * 5..NWORDS_FIELD * 6]);
}

pub fn fp2_encode(x: &f2elm_t, enc: &mut [u8]) {
    let mut t: f2elm_t = f2elm;
    fpx::from_fp2mont(x, &mut t);

    let a = F2ELM::from(t).0;

    for i in 0..(crate::FP2_ENCODED_BYTES / 2) {
        enc[i] = a[i];
        enc[i + crate::FP2_ENCODED_BYTES / 2] = a[i + crate::MAXBITS_FIELD / 8];
    }
}

pub fn fp2_decode(enc: &[u8], x: &mut f2elm_t) {
    let mut a = F2ELM::from(*x);

    for i in 0..2 * (crate::MAXBITS_FIELD / 8) {
        a.0[i] = 0;
    }

    for i in 0..crate::FP2_ENCODED_BYTES / 2 {
        a.0[i] = enc[i];
        a.0[i + crate::MAXBITS_FIELD / 8] = enc[i + crate::FP2_ENCODED_BYTES / 2];
    }

    let a = f2elm_t::from(a);

    fpx::to_fp2mont(&a, x);
}

pub fn random_mod_order_A(random_digits: &mut [u8]) {
    let nbytes = (crate::OALICE_BITS + 7) / 8;

    random_digits[..crate::MAXWORDS_ORDER].zeroize();
    getrandom::getrandom(&mut random_digits[..nbytes]).unwrap();
    random_digits[nbytes - 1] = (random_digits[nbytes - 1] as i32 & crate::MASK_ALICE as i32) as u8;
}

pub fn random_mod_order_B(random_digits: &mut [u8]) {
    let nbytes = (crate::OBOB_BITS - 1 + 7) / 8;

    random_digits[..crate::MAXWORDS_ORDER].zeroize();
    getrandom::getrandom(&mut random_digits[..nbytes]).unwrap();
    random_digits[nbytes - 1] = (random_digits[nbytes - 1] as i32 & crate::MASK_ALICE as i32) as u8;
}

pub fn ephemeral_key_generation_A(private_key_a: &[u8], public_key_a: &mut [u8]) -> i32 {
    let mut r = [PointProj::new()];
    let mut phiP = PointProj::new();
    let mut phiQ = PointProj::new();
    let mut phiR = PointProj::new();
    let mut pts = [PointProj::new(); 8];
    let mut XPA = f2elm;
    let mut XQA = f2elm;
    let mut XRA = f2elm;
    let mut coeff: [f2elm_t; 3] = [f2elm; 3];
    let mut a24plus = f2elm;
    let mut c24 = f2elm;
    let mut a = f2elm;
    let mut pts_index: [libc::c_uint; 8] = [0; 8];
    let mut npts: libc::c_uint = 0 as libc::c_int as libc::c_uint;
    let mut ii: libc::c_uint = 0 as libc::c_int as libc::c_uint;

    init_basis(&crate::A_gen, &mut XPA, &mut XQA, &mut XRA);
    init_basis(&crate::B_gen, &mut phiP.x, &mut phiQ.x, &mut phiR.x);

    phiP.z[0] = crate::Montgomery_one;
    phiQ.z[0] = crate::Montgomery_one;
    phiR.z[0] = crate::Montgomery_one;
    a24plus[0] = crate::Montgomery_one;

    fpx::fp2addx(a24plus.as_ptr(), a24plus.as_ptr(), &mut a24plus);
    fpx::fp2add(&a24plus, &a24plus, &mut c24);
    fpx::fp2add(&a24plus, &c24, &mut a);
    fpx::fp2add(&c24, &c24, &mut a24plus);

    ec_isogeny::ladder3pt(
        &XPA,
        &XQA,
        &XRA,
        private_key_a.as_ptr() as *const u64,
        0,
        &mut r[0],
        &a,
    );

    let mut index: u32 = 0;

    for row in 1..crate::MAX_Alice as u32 {
        while index < (crate::MAX_Alice as u32).wrapping_sub(row) {
            fpx::fp2copy(&r[0].x, &mut pts[npts as usize].x);
            fpx::fp2copy(&r[0].z, &mut pts[npts as usize].z);

            let fresh2 = npts;
            npts = npts.wrapping_add(1);
            pts_index[fresh2 as usize] = index;
            let fresh3 = ii;
            ii = ii.wrapping_add(1);
            let m = crate::strat_Alice[fresh3 as usize];

            ec_isogeny::xDBLe(
                r.as_ptr(),
                &mut r[0],
                &a24plus,
                &c24,
                (2 as u32).wrapping_mul(m) as usize,
            );

            index = index.wrapping_add(m);
        }

        ec_isogeny::get_4_isog(&r[0], &mut a24plus, &mut c24, &mut coeff);

        for i in 0..npts {
            ec_isogeny::eval_4_isog(&mut pts[i as usize], &mut coeff);
        }

        ec_isogeny::eval_4_isog(&mut phiP, &mut coeff);
        ec_isogeny::eval_4_isog(&mut phiQ, &mut coeff);
        ec_isogeny::eval_4_isog(&mut phiR, &mut coeff);
        fpx::fp2copy(&pts[npts.wrapping_sub(1) as usize].x, &mut r[0].x);
        fpx::fp2copy(&pts[npts.wrapping_sub(1) as usize].z, &mut r[0].z);

        index = pts_index[npts.wrapping_sub(1) as usize];
        npts = npts.wrapping_sub(1);
    }

    ec_isogeny::get_4_isog(&r[0], &mut a24plus, &mut c24, &mut coeff);
    ec_isogeny::eval_4_isog(&mut phiP, &mut coeff);
    ec_isogeny::eval_4_isog(&mut phiQ, &mut coeff);
    ec_isogeny::eval_4_isog(&mut phiR, &mut coeff);
    ec_isogeny::inv_3_way(&mut phiP.z, &mut phiQ.z, &mut phiR.z);
    fpx::fp2mul_montx(phiP.x.as_ptr(), phiP.z.as_ptr(), &mut phiP.x);
    fpx::fp2mul_montx(phiQ.x.as_ptr(), phiQ.z.as_ptr(), &mut phiQ.x);
    fpx::fp2mul_montx(phiR.x.as_ptr(), phiR.z.as_ptr(), &mut phiR.x);

    // Format public key
    fp2_encode(&phiP.x, public_key_a);
    fp2_encode(&phiQ.x, &mut public_key_a[crate::FP2_ENCODED_BYTES..]);
    fp2_encode(&phiR.x, &mut public_key_a[2 * crate::FP2_ENCODED_BYTES..]);

    return 0;
}

#[cfg(test)]
pub mod tests {
    use crate::ec_isogeny::tests::*;
    use crate::fp_generic::tests::*;
    use crate::fpx::tests::*;
    use crate::util::fillbytes_u64;
    use crate::*;

    #[test]
    fn test_sidh() {
        let a = [
            251, 164, 62, 62, 97, 116, 208, 159, 111, 107, 17, 56, 93, 111, 187, 150, 211, 43, 202,
            244, 134, 226, 173, 216, 24, 21, 199, 28, 22, 117, 38, 114, 43, 29, 132, 68, 209, 247,
            166, 111, 178, 161, 56, 150, 182, 177, 244, 50, 102, 23, 165, 252, 125, 21, 189, 140,
            93, 76, 245, 238, 246, 49, 224, 160, 110, 202, 147, 253, 71, 248, 160, 208, 89, 3, 174,
            6, 95, 136, 223, 87, 213, 62, 244, 56, 135, 255, 182, 38, 155, 233, 97, 171, 179, 85,
            250, 148, 25, 33, 19, 66, 21, 251, 26, 87, 27, 35, 0, 12, 180, 205, 39, 226, 9, 199,
            11, 51, 213, 91, 196, 86, 19, 205, 8, 25, 154, 39, 66, 217, 148, 139, 17, 113, 239,
            199, 241, 79, 17, 62, 144, 113, 124, 167, 216, 70, 216, 8, 170, 122, 173, 4, 80, 207,
            182, 114, 157, 158, 96, 193, 162, 48, 215, 32, 68, 165, 11, 33, 96, 221, 35, 58, 177,
            140, 219, 160, 116, 203, 183, 4, 195, 42, 220, 174, 174, 83, 191, 135, 212, 117, 202,
            37, 33, 18, 197, 92, 163, 162, 73, 87, 195, 178, 61, 72, 184, 245, 115, 144, 35, 69,
            10, 169, 63, 151, 38, 72, 105, 179, 227, 134, 254, 197, 204, 249, 175, 134, 105, 83,
            233, 51, 130, 175, 3, 150, 31, 59, 219, 219, 68, 91, 57, 94, 175, 229, 51, 84, 123,
            249, 35, 187, 8, 169, 97, 65, 116, 109, 139, 15, 208, 8, 8, 69, 71, 161, 211, 177, 6,
            120, 84, 32, 102, 219, 63, 235, 110, 0, 171, 144, 156, 181, 138, 213, 148, 191, 23,
            116, 157, 108, 201, 54, 133, 220, 204, 22, 224, 154, 21, 127, 252, 117, 72, 126, 25,
            38, 177, 222, 18, 15, 70, 115, 230, 136, 238, 138, 182, 154, 72, 71, 62, 96, 15, 236,
            156, 86, 99, 174, 55, 93, 92, 121, 123, 118, 191, 251, 90, 246, 240, 244, 111, 141,
            213, 177, 179, 215, 23, 248, 111, 207, 194, 239, 179, 88, 34, 73, 118, 161, 95, 96,
            171, 141, 211, 255, 67, 6, 123, 32, 196, 124, 118, 114, 235, 23, 168, 186, 46, 130,
            143, 52, 69, 125, 131, 69, 232, 29, 248, 233, 143, 108, 85, 2, 136, 202, 123, 201, 42,
            238, 44, 37, 111, 124, 99, 15, 27, 109, 183, 5, 174, 169, 184, 22, 67, 114, 93, 0, 17,
            184, 240, 27, 164, 99, 27, 6, 11, 59, 35, 208, 50, 149, 130, 153, 63, 20, 25, 142, 154,
            99, 224, 244, 169, 25, 173, 33, 185, 140, 40, 169, 103, 73, 193, 37, 251, 13, 199, 192,
            253, 51, 115, 192, 54, 227, 238, 64, 248, 148, 106, 100, 6, 18, 235, 96, 222, 186, 186,
            130, 47, 143, 78, 13, 233, 47, 116, 253, 124, 114, 71, 50, 251, 253, 32, 64, 182, 172,
            169, 234, 217, 143, 43, 60, 88, 81, 181, 60, 216, 75, 13, 227, 254, 27, 202, 182, 193,
            255, 197, 138, 2, 35, 63, 168, 229, 57, 151, 66, 92, 3, 10, 253, 67, 137, 105, 56, 111,
            149, 21, 5, 178, 98, 250, 157, 0, 184, 243, 30, 158, 149, 46, 174, 185, 75, 77, 249,
            104, 194, 74, 4, 118, 123, 32, 102, 36, 181, 110, 219, 110, 35, 38, 155, 221, 161, 174,
            210, 79, 181, 89, 49, 121, 73, 251, 227, 83, 87, 249, 42, 129, 90, 128, 16, 120, 137,
            84, 96, 70, 30, 109, 240, 28, 246, 145, 231, 31, 3, 38, 61, 11, 88, 43, 225, 51, 24,
            56, 104, 151, 215, 30, 57, 193, 22, 254, 252, 96, 78, 183, 227, 184, 12, 15, 0, 95,
            167, 67, 183, 137, 169, 218, 235, 218, 235, 45, 119, 38, 153, 110, 196, 177, 11, 125,
            50, 120, 43, 28, 83, 143, 151, 25, 138, 129, 52, 92, 151, 31, 185, 223, 174, 88, 156,
            95, 241, 38, 193, 207, 85, 237, 11, 139, 233, 233, 97, 155, 232, 21, 87, 212, 12, 142,
            213, 91, 219, 145, 188, 243, 253, 27, 224, 124, 227, 85, 71, 250, 210, 144, 105, 172,
            32, 71, 91, 17, 17, 18, 197, 31, 104, 64, 116, 83, 205, 246, 209, 161, 192, 51, 26,
            178, 148, 180, 103, 186, 217, 80, 15, 152, 92, 136, 236, 189, 153, 215, 193, 125, 131,
            229, 6, 65, 218, 19, 245, 209, 218, 156, 66, 130, 60, 147, 255, 152, 145, 90, 74, 68,
            63, 2, 53, 62, 155, 242, 207, 9, 92, 27, 120, 66, 216, 172, 167, 232, 88, 83, 243, 30,
            190, 116, 72, 85, 255, 93, 23, 22, 109, 234, 128, 221, 34, 9, 31, 254, 185, 18, 28,
            129, 172, 240, 97, 146, 188, 227, 61, 0, 184, 5, 247, 150, 211, 31, 240, 69, 110, 36,
            222, 193, 219, 89, 231, 21, 77, 187, 175, 237, 30, 17, 180, 152, 180, 89, 28, 238, 148,
            12, 72, 25, 195, 16, 227, 16, 216, 205, 192, 94, 229, 57, 221, 11, 163, 174, 232, 14,
            210, 116, 60, 106, 161, 209, 117, 167, 215, 211, 162, 139, 64, 14, 214, 146, 239, 33,
            141, 254, 78, 229, 51, 55, 114, 250, 129, 91, 21, 190, 18, 234, 230, 17, 2, 254, 120,
            145, 80, 216, 129, 48, 241, 38, 149, 225, 149, 131, 243, 33, 99, 27, 246, 10, 132, 228,
            248, 154, 214, 42, 126, 79, 162, 166, 206, 15, 210, 138, 52, 56, 148, 10, 64, 116, 131,
            34, 11, 170, 143, 220, 58, 182, 169, 225, 2, 225, 233, 80, 196, 181, 185, 161, 249, 85,
            24, 24, 32, 239, 29, 86, 52, 73, 110, 104, 225, 165, 191, 45, 204, 119, 15, 38, 72, 17,
            204, 70, 254, 176, 102, 214, 241, 22, 15, 10, 237, 53, 160, 230, 239, 140, 136, 201,
            123, 238, 208, 31, 198, 94, 83, 157, 23, 130, 157, 91, 81, 187, 230, 209, 137,
        ];
        let mut b = [0u8; 564];

        sidh::ephemeral_key_generation_A(&a, &mut b);

        assert_eq!(
            b.as_ref(),
            [
                55, 187, 16, 127, 0, 231, 216, 34, 158, 62, 136, 249, 145, 43, 78, 168, 90, 205,
                115, 62, 173, 218, 148, 50, 193, 60, 167, 218, 238, 228, 214, 121, 114, 137, 229,
                201, 29, 143, 188, 14, 106, 201, 214, 122, 93, 9, 212, 39, 188, 193, 107, 191, 191,
                218, 164, 170, 62, 138, 68, 6, 208, 202, 207, 235, 221, 14, 228, 225, 66, 87, 40,
                175, 202, 134, 41, 121, 16, 82, 221, 17, 182, 176, 180, 248, 114, 150, 93, 70, 38,
                23, 102, 154, 99, 39, 2, 37, 207, 45, 52, 46, 191, 33, 49, 238, 63, 245, 232, 238,
                217, 231, 84, 10, 251, 48, 219, 110, 238, 184, 201, 223, 212, 58, 80, 27, 65, 163,
                152, 135, 6, 134, 250, 177, 148, 66, 65, 190, 208, 109, 160, 133, 32, 127, 128,
                153, 59, 82, 196, 130, 115, 55, 66, 32, 54, 141, 222, 179, 245, 207, 9, 141, 134,
                14, 49, 122, 141, 102, 35, 53, 142, 210, 221, 17, 186, 193, 171, 178, 206, 226,
                243, 36, 199, 110, 76, 141, 83, 127, 69, 15, 229, 97, 245, 60, 194, 212, 115, 33,
                160, 81, 168, 36, 198, 97, 119, 134, 228, 241, 15, 251, 110, 89, 54, 51, 180, 43,
                6, 72, 50, 89, 137, 189, 21, 9, 82, 197, 209, 129, 128, 6, 149, 200, 136, 148, 68,
                248, 250, 17, 193, 74, 12, 49, 225, 82, 102, 94, 22, 124, 73, 6, 129, 18, 104, 246,
                151, 217, 196, 251, 224, 226, 12, 245, 217, 127, 160, 191, 239, 88, 177, 109, 105,
                159, 226, 192, 45, 157, 240, 155, 189, 157, 146, 153, 169, 54, 69, 199, 12, 175,
                255, 74, 3, 181, 234, 206, 47, 85, 122, 67, 48, 241, 161, 166, 93, 23, 62, 57, 104,
                59, 198, 126, 251, 144, 33, 112, 234, 242, 246, 252, 164, 165, 138, 43, 221, 110,
                30, 163, 64, 88, 70, 72, 233, 52, 131, 140, 244, 169, 222, 146, 63, 253, 61, 131,
                44, 204, 107, 165, 23, 240, 240, 205, 230, 167, 201, 30, 151, 200, 125, 185, 41,
                198, 145, 29, 251, 226, 101, 209, 233, 22, 43, 33, 64, 94, 25, 194, 164, 171, 79,
                11, 21, 215, 202, 206, 66, 105, 181, 58, 145, 241, 166, 209, 209, 91, 27, 250, 110,
                15, 239, 105, 75, 208, 146, 235, 109, 49, 230, 229, 251, 226, 4, 160, 181, 177, 79,
                216, 2, 132, 29, 73, 31, 171, 186, 169, 98, 122, 66, 216, 46, 246, 241, 241, 54,
                108, 44, 243, 139, 36, 243, 185, 157, 199, 193, 128, 62, 219, 174, 93, 66, 68, 73,
                118, 134, 130, 215, 117, 220, 33, 199, 36, 13, 227, 140, 233, 132, 28, 212, 75,
                168, 17, 53, 141, 70, 33, 57, 238, 243, 141, 51, 205, 57, 189, 12, 189, 196, 194,
                168, 149, 245, 221, 118, 146, 179, 175, 154, 231, 57, 218, 181, 133, 103, 245, 199,
                248, 222, 32, 92, 230, 230, 244, 187, 164, 147, 173, 254, 251, 153, 72, 36, 138,
                142, 49, 247, 66, 38, 2, 239, 104, 108, 251, 248, 3, 135, 228, 131, 203, 101, 248,
                51, 224, 107, 167, 114, 136, 52, 232, 22, 133, 211, 42, 116, 53, 16, 67, 50, 183,
                47, 250, 196, 4, 198, 206, 120, 51, 140, 178, 195, 55
            ]
            .as_ref()
        );

        for _ in 0..500 {
            unsafe {
                let mut x = [[0u64; 12]; 2];
                let mut enc = [0u8; 10000];
                let mut end = [0u8; 10000];

                fillbytes_u64(&mut [&mut x[0]]);
                fillbytes_u64(&mut [&mut x[1]]);

                sidh::fp2_encode(&x, &mut enc);
                fp2_encode(x.as_ptr(), end.as_mut_ptr());

                assert_eq!(enc.as_ref(), end.as_ref());

                let mut xa = [[0u64; 12]; 2];
                let mut xb = [[0u64; 12]; 2];

                sidh::fp2_decode(&enc, &mut xa);
                fp2_decode(end.as_ptr(), xb.as_mut_ptr());

                //assert_eq!(xa, xb);
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn init_basis(
        mut gen: *mut digit_t,
        mut XP: *mut felm_t,
        mut XQ: *mut felm_t,
        mut XR: *mut felm_t,
    ) {
        // Initialization of basis points
        fpcopy(
            gen as *const digit_t,
            (*XP.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpcopy(
            gen.offset(12 as libc::c_int as isize) as *const digit_t,
            (*XP.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpcopy(
            gen.offset((2 as libc::c_int * 12 as libc::c_int) as isize) as *const digit_t,
            (*XQ.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpcopy(
            gen.offset((3 as libc::c_int * 12 as libc::c_int) as isize) as *const digit_t,
            (*XQ.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpcopy(
            gen.offset((4 as libc::c_int * 12 as libc::c_int) as isize) as *const digit_t,
            (*XR.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpcopy(
            gen.offset((5 as libc::c_int * 12 as libc::c_int) as isize) as *const digit_t,
            (*XR.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn EphemeralKeyGeneration_A(
        mut PrivateKeyA: *const libc::c_uchar,
        mut PublicKeyA: *mut libc::c_uchar,
    ) -> libc::c_int {
        // Alice's ephemeral public key generation
        // Input:  a private key PrivateKeyA in the range [0, 2^eA - 1].
        // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
        let mut R: point_proj_t = [point_proj { X: f2elm, Z: f2elm }; 1];
        let mut phiP: point_proj_t = [{
            let mut init = point_proj {
                X: [
                    [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0; 12],
                ],
                Z: f2elm,
            };
            init
        }];
        let mut phiQ: point_proj_t = [{
            let mut init = point_proj {
                X: [
                    [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0; 12],
                ],
                Z: f2elm,
            };
            init
        }];
        let mut phiR: point_proj_t = [{
            let mut init = point_proj {
                X: [
                    [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0; 12],
                ],
                Z: f2elm,
            };
            init
        }];
        let mut pts: [point_proj_t; 8] = [[point_proj { X: f2elm, Z: f2elm }; 1]; 8];
        let mut XPA: f2elm_t = f2elm;
        let mut XQA: f2elm_t = f2elm;
        let mut XRA: f2elm_t = f2elm;
        let mut coeff: [f2elm_t; 3] = [f2elm; 3];
        let mut A24plus: f2elm_t = [
            [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0; 12],
        ];
        let mut C24: f2elm_t = [
            [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0; 12],
        ];
        let mut A: f2elm_t = [
            [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0; 12],
        ];
        let mut i: libc::c_uint = 0;
        let mut row: libc::c_uint = 0;
        let mut m: libc::c_uint = 0;
        let mut index: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        let mut pts_index: [libc::c_uint; 8] = [0; 8];
        let mut npts: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        let mut ii: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        // Initialize basis points
        init_basis(
            A_gen.as_ptr() as *mut digit_t,
            XPA.as_mut_ptr(),
            XQA.as_mut_ptr(),
            XRA.as_mut_ptr(),
        );
        init_basis(
            B_gen.as_ptr() as *mut digit_t,
            (*phiP.as_mut_ptr()).X.as_mut_ptr(),
            (*phiQ.as_mut_ptr()).X.as_mut_ptr(),
            (*phiR.as_mut_ptr()).X.as_mut_ptr(),
        );
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            (*phiP.as_mut_ptr()).Z[0 as libc::c_int as usize].as_mut_ptr(),
        );
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            (*phiQ.as_mut_ptr()).Z[0 as libc::c_int as usize].as_mut_ptr(),
        );
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            (*phiR.as_mut_ptr()).Z[0 as libc::c_int as usize].as_mut_ptr(),
        );
        // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            A24plus[0 as libc::c_int as usize].as_mut_ptr(),
        );
        fp2add(
            A24plus.as_mut_ptr() as *const felm_t,
            A24plus.as_mut_ptr() as *const felm_t,
            A24plus.as_mut_ptr(),
        );
        fp2add(
            A24plus.as_mut_ptr() as *const felm_t,
            A24plus.as_mut_ptr() as *const felm_t,
            C24.as_mut_ptr(),
        );
        fp2add(
            A24plus.as_mut_ptr() as *const felm_t,
            C24.as_mut_ptr() as *const felm_t,
            A.as_mut_ptr(),
        );
        fp2add(
            C24.as_mut_ptr() as *const felm_t,
            C24.as_mut_ptr() as *const felm_t,
            A24plus.as_mut_ptr(),
        );
        // Retrieve kernel point
        LADDER3PT(
            XPA.as_mut_ptr() as *const felm_t,
            XQA.as_mut_ptr() as *const felm_t,
            XRA.as_mut_ptr() as *const felm_t,
            PrivateKeyA as *mut digit_t,
            0 as libc::c_int as libc::c_uint,
            R.as_mut_ptr(),
            A.as_mut_ptr() as *const felm_t,
        );
        // Traverse tree
        index = 0 as libc::c_int as libc::c_uint;
        row = 1 as libc::c_int as libc::c_uint;
        while row < 186 as libc::c_int as libc::c_uint {
            while index < (186 as libc::c_int as libc::c_uint).wrapping_sub(row) {
                fp2copy(
                    (*R.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
                    (*pts[npts as usize].as_mut_ptr()).X.as_mut_ptr(),
                );
                fp2copy(
                    (*R.as_mut_ptr()).Z.as_mut_ptr() as *const felm_t,
                    (*pts[npts as usize].as_mut_ptr()).Z.as_mut_ptr(),
                );
                let fresh2 = npts;
                npts = npts.wrapping_add(1);
                pts_index[fresh2 as usize] = index;
                let fresh3 = ii;
                ii = ii.wrapping_add(1);
                m = strat_Alice[fresh3 as usize];
                xDBLe(
                    R.as_mut_ptr() as *const point_proj,
                    R.as_mut_ptr(),
                    A24plus.as_mut_ptr() as *const felm_t,
                    C24.as_mut_ptr() as *const felm_t,
                    (2 as libc::c_int as libc::c_uint).wrapping_mul(m) as libc::c_int,
                );
                index = index.wrapping_add(m)
            }
            get_4_isog(
                R.as_mut_ptr() as *const point_proj,
                A24plus.as_mut_ptr(),
                C24.as_mut_ptr(),
                coeff.as_mut_ptr(),
            );
            i = 0 as libc::c_int as libc::c_uint;
            while i < npts {
                eval_4_isog(pts[i as usize].as_mut_ptr(), coeff.as_mut_ptr());
                i = i.wrapping_add(1)
            }
            eval_4_isog(phiP.as_mut_ptr(), coeff.as_mut_ptr());
            eval_4_isog(phiQ.as_mut_ptr(), coeff.as_mut_ptr());
            eval_4_isog(phiR.as_mut_ptr(), coeff.as_mut_ptr());
            fp2copy(
                (*pts[npts.wrapping_sub(1 as libc::c_int as libc::c_uint) as usize].as_mut_ptr())
                    .X
                    .as_mut_ptr() as *const felm_t,
                (*R.as_mut_ptr()).X.as_mut_ptr(),
            );
            fp2copy(
                (*pts[npts.wrapping_sub(1 as libc::c_int as libc::c_uint) as usize].as_mut_ptr())
                    .Z
                    .as_mut_ptr() as *const felm_t,
                (*R.as_mut_ptr()).Z.as_mut_ptr(),
            );
            index = pts_index[npts.wrapping_sub(1 as libc::c_int as libc::c_uint) as usize];
            npts = npts.wrapping_sub(1 as libc::c_int as libc::c_uint);
            row = row.wrapping_add(1)
        }
        get_4_isog(
            R.as_mut_ptr() as *const point_proj,
            A24plus.as_mut_ptr(),
            C24.as_mut_ptr(),
            coeff.as_mut_ptr(),
        );
        eval_4_isog(phiP.as_mut_ptr(), coeff.as_mut_ptr());
        eval_4_isog(phiQ.as_mut_ptr(), coeff.as_mut_ptr());
        eval_4_isog(phiR.as_mut_ptr(), coeff.as_mut_ptr());
        inv_3_way(
            (*phiP.as_mut_ptr()).Z.as_mut_ptr(),
            (*phiQ.as_mut_ptr()).Z.as_mut_ptr(),
            (*phiR.as_mut_ptr()).Z.as_mut_ptr(),
        );
        fp2mul_mont(
            (*phiP.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
            (*phiP.as_mut_ptr()).Z.as_mut_ptr() as *const felm_t,
            (*phiP.as_mut_ptr()).X.as_mut_ptr(),
        );
        fp2mul_mont(
            (*phiQ.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
            (*phiQ.as_mut_ptr()).Z.as_mut_ptr() as *const felm_t,
            (*phiQ.as_mut_ptr()).X.as_mut_ptr(),
        );
        fp2mul_mont(
            (*phiR.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
            (*phiR.as_mut_ptr()).Z.as_mut_ptr() as *const felm_t,
            (*phiR.as_mut_ptr()).X.as_mut_ptr(),
        );
        // Format public key
        fp2_encode(
            (*phiP.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
            PublicKeyA,
        );
        fp2_encode(
            (*phiQ.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
            PublicKeyA.offset(
                (2 as libc::c_int * ((751 as libc::c_int + 7 as libc::c_int) / 8 as libc::c_int))
                    as isize,
            ),
        );
        fp2_encode(
            (*phiR.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
            PublicKeyA.offset(
                (2 as libc::c_int
                    * 2 as libc::c_int
                    * ((751 as libc::c_int + 7 as libc::c_int) / 8 as libc::c_int))
                    as isize,
            ),
        );
        return 0 as libc::c_int;
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2_encode(mut x: *const felm_t, mut enc: *mut libc::c_uchar) {
        // Conversion of GF(p^2) element from Montgomery to standard representation, and encoding by removing leading 0 bytes
        let mut i: libc::c_uint = 0;
        let mut t: f2elm_t = f2elm;
        from_fp2mont(x, t.as_mut_ptr());
        i = 0 as libc::c_int as libc::c_uint;
        while i
            < (2 as libc::c_int * ((751 as libc::c_int + 7 as libc::c_int) / 8 as libc::c_int)
                / 2 as libc::c_int) as libc::c_uint
        {
            *enc.offset(i as isize) = *(t.as_mut_ptr() as *mut libc::c_uchar).offset(i as isize);
            *enc.offset(i.wrapping_add(
                (2 as libc::c_int * ((751 as libc::c_int + 7 as libc::c_int) / 8 as libc::c_int)
                    / 2 as libc::c_int) as libc::c_uint,
            ) as isize) = *(t.as_mut_ptr() as *mut libc::c_uchar).offset(
                i.wrapping_add((768 as libc::c_int / 8 as libc::c_int) as libc::c_uint) as isize,
            );
            i = i.wrapping_add(1)
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2_decode(mut enc: *const libc::c_uchar, mut x: *mut felm_t) {
        // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation
        let mut i: libc::c_uint = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < (2 as libc::c_int * (768 as libc::c_int / 8 as libc::c_int)) as libc::c_uint {
            *(x as *mut libc::c_uchar).offset(i as isize) = 0 as libc::c_int as libc::c_uchar;
            i = i.wrapping_add(1)
        }
        i = 0 as libc::c_int as libc::c_uint;
        while i
            < (2 as libc::c_int * ((751 as libc::c_int + 7 as libc::c_int) / 8 as libc::c_int)
                / 2 as libc::c_int) as libc::c_uint
        {
            *(x as *mut libc::c_uchar).offset(i as isize) = *enc.offset(i as isize);
            *(x as *mut libc::c_uchar).offset(
                i.wrapping_add((768 as libc::c_int / 8 as libc::c_int) as libc::c_uint) as isize,
            ) = *enc.offset(i.wrapping_add(
                (2 as libc::c_int * ((751 as libc::c_int + 7 as libc::c_int) / 8 as libc::c_int)
                    / 2 as libc::c_int) as libc::c_uint,
            ) as isize);
            i = i.wrapping_add(1)
        }
        to_fp2mont(x as *const felm_t, x);
    }
}
