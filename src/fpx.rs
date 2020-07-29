use crate::util::{addc, shiftl, subc};
use crate::{f2elm_t, felm_t};

pub fn to_mont(a: &felm_t, mc: &mut felm_t) {
    fpmul_mont(a, &crate::Montgomery_R2, mc);
}

pub fn from_mont(ma: &felm_t, c: &mut felm_t) {
    let mut one = [0u64; crate::NWORDS_FIELD];
    one[0] = 1;

    fpmul_mont(ma, &one, c);
    crate::fp_generic::fpcorrection(c);
}

pub fn copy_words(a: &[u64], c: &mut [u64], nwords: usize) {
    for i in 0..nwords {
        c[i] = a[i];
    }
}

pub fn fpmul_mont(ma: &felm_t, mb: &felm_t, mc: &mut felm_t) {
    let mut temp = [0u64; 24];

    crate::fp_generic::mp_mul(ma, mb, &mut temp, 12);
    crate::fp_generic::rdc_mont(&temp, mc);
}

pub fn fpsqr_mont(ma: &felm_t, mc: &mut felm_t) {
    let mut temp = [0u64; 24];

    crate::fp_generic::mp_mul(ma, ma, &mut temp, 12);
    crate::fp_generic::rdc_mont(&temp, mc);
}

pub fn fpinv_mont(a: &mut felm_t) {
    let mut tt = *a;
    let mut temp = tt;

    fpinv_chain_mont(a);
    fpsqr_mont(&temp, &mut tt);
    temp = tt;
    fpsqr_mont(&temp, &mut tt);
    temp = *a;
    fpmul_mont(&temp, &tt, a);
    todo!();
}

pub fn fp2copy(a: &f2elm_t, c: &mut f2elm_t) {
    c[0] = a[0];
    c[1] = a[1];
}

pub fn fp2zero(a: &mut f2elm_t) {
    a[0].iter_mut().for_each(|x| *x = 0);
    a[1].iter_mut().for_each(|x| *x = 0);
}

pub fn fp2neg(a: &mut f2elm_t) {
    crate::fp_generic::fpneg(&mut a[0]);
    crate::fp_generic::fpneg(&mut a[1]);
}

pub fn fp2add(a: &f2elm_t, b: &f2elm_t, c: &mut f2elm_t) {
    crate::fp_generic::fpadd(&a[0], &b[0], &mut c[0]);
    crate::fp_generic::fpadd(&a[1], &b[1], &mut c[1]);
}

pub fn fp2sub(a: &f2elm_t, b: &f2elm_t, c: &mut f2elm_t) {
    crate::fp_generic::fpsub(&a[0], &b[0], &mut c[0]);
    crate::fp_generic::fpsub(&a[1], &b[1], &mut c[1]);
}

pub fn fp2div2(a: &f2elm_t, c: &mut f2elm_t) {
    crate::fp_generic::fpdiv2(&a[0], &mut c[0]);
}

pub fn fp2correction(a: &mut f2elm_t) {
    crate::fp_generic::fpcorrection(&mut a[0]);
    crate::fp_generic::fpcorrection(&mut a[1]);
}

pub fn mp_addfast(a: &[u64], b: &[u64], c: &mut [u64]) {
    mp_add(a, b, c, 12);
}

pub fn fp2sqr_mont(a: &f2elm_t, c: &mut f2elm_t) {
    let mut t1 = [0u64; 12];
    let mut t2 = [0u64; 12];
    let mut t3 = [0u64; 12];

    mp_addfast(&a[0], &a[1], &mut t1);
    crate::fp_generic::fpsub(&a[0], &a[1], &mut t2);
    mp_addfast(&a[0], &a[0], &mut t3);
    fpmul_mont(&t1, &t2, &mut c[0]);
    fpmul_mont(&t3, &a[1], &mut c[1]);
}

pub fn mp_sub(a: &[u64], b: &[u64], c: &mut [u64], nwords: usize) -> u32 {
    let mut borrow = 0;

    for i in 0..nwords {
        subc(borrow, a[i], b[i], &mut borrow, &mut c[i]);
    }

    borrow
}

pub fn mp_subfast(a: &[u64], b: &[u64], c: &mut [u64]) -> u64 {
    0 - mp_sub(a, b, c, 24) as u64
}

pub fn mp_dblsubfast(a: &[u64], b: &[u64], c: &mut [u64]) {
    let mut temp = c.to_vec();

    mp_sub(c, a, &mut temp, 24);
    c.clone_from_slice(&temp);

    mp_sub(c, b, &mut temp, 24);
    c.clone_from_slice(&temp);
}

pub fn fp2mul_mont(a: &f2elm_t, b: &f2elm_t, c: &mut f2elm_t) {
    let mut t1 = [0u64; 12];
    let mut t2 = [0u64; 12];
    let mut tt1 = [0u64; 24];
    let mut tt2 = [0u64; 24];
    let mut tt3 = [0u64; 24];

    mp_addfast(&a[0], &a[1], &mut t1);
    mp_addfast(&b[0], &b[1], &mut t2);
    crate::fp_generic::mp_mul(&a[0], &b[0], &mut tt1, 12);
    crate::fp_generic::mp_mul(&a[1], &b[1], &mut tt2, 12);
    crate::fp_generic::mp_mul(&t1, &t2, &mut tt3, 12);
    mp_dblsubfast(&tt1, &tt2, &mut tt3);

    let mut temp = tt1;
    let mask = mp_subfast(&tt1, &tt2, &mut temp);
    tt1 = temp;

    for i in 0..12 as usize {
        t1[i] = crate::p751[i] & mask;
    }

    crate::fp_generic::rdc_mont(&tt3, &mut c[1]);
    let mut temp = tt1;
    mp_addfast(&tt1[12..], &t1, &mut temp[12..]);
    tt1 = temp;
    crate::fp_generic::rdc_mont(&tt1, &mut c[0]);
}

pub fn fpinv_chain_mont(a: &mut felm_t) {
    let mut t = [[0u64; 12]; 27];
    let mut tt = [0u64; 12];

    fpsqr_mont(a, &mut tt);
    fpmul_mont(a, &mut tt, &mut t[0]);

    let temp = t[0];
    fpmul_mont(&temp, &tt, &mut t[1]);
    let temp = t[0];
    fpmul_mont(&temp, &tt, &mut t[2]);
    let temp = t[0];
    fpmul_mont(&temp, &tt, &mut t[3]);
    let temp = t[0];
    fpmul_mont(&temp, &tt, &mut t[3]);

    for i in 3..9 {
        let temp = t[i];
        fpmul_mont(&temp, &tt, &mut t[i + 1]);
    }

    let temp = t[9];
    fpmul_mont(&temp, &tt, &mut t[9]);

    for i in 9..21 {
        let temp = t[i];
        fpmul_mont(&temp, &tt, &mut t[i + 1]);
    }

    let temp = t[21];
    fpmul_mont(&temp, &tt, &mut t[21]);

    for i in 21..25 {
        let temp = t[i];
        fpmul_mont(&temp, &tt, &mut t[i + 1]);
    }

    let temp = t[25];
    fpmul_mont(&temp, &tt, &mut t[25]);
    let temp = t[25];
    fpmul_mont(&temp, &tt, &mut t[26]);

    tt = *a;

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[20], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[24], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[11], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[8], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[2], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[23], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[2], &temp, &mut tt);

    for _ in 0..9 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[2], &temp, &mut tt);

    for _ in 0..10 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[15], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[13], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[26], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[20], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[11], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[10], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[14], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[4], &temp, &mut tt);

    for _ in 0..10 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[18], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[1], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[22], &temp, &mut tt);

    for _ in 0..10 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[6], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[24], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[9], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[18], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[17], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(a, &temp, &mut tt);

    for _ in 0..10 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[16], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[7], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[0], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[12], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[19], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[22], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[25], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[2], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[10], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[22], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[18], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[4], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[14], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[13], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[5], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[23], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[21], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[2], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[23], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[12], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[9], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[3], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[13], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[17], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[26], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[5], &temp, &mut tt);

    for _ in 0..8 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[8], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[2], &temp, &mut tt);

    for _ in 0..6 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[11], &temp, &mut tt);

    for _ in 0..7 {
        let temp = tt;
        fpsqr_mont(&temp, &mut tt);
    }

    let temp = tt;
    fpmul_mont(&t[20], &temp, &mut tt);

    for _ in 0..61 {
        for _ in 0..6 {
            let temp = tt;
            fpsqr_mont(&temp, &mut tt);
        }
        let temp = tt;
        fpmul_mont(&t[26], &temp, &mut tt);
    }

    *a = tt;
}

pub fn fp2inv_mont(a: &mut f2elm_t) {
    let mut t1 = [[0u64; 12]; 2];

    fpsqr_mont(&a[0], &mut t1[0]);
    fpsqr_mont(&a[1], &mut t1[1]);

    let temp = t1;
    crate::fp_generic::fpadd(&temp[0], &temp[1], &mut t1[0]);
    fpinv_mont(&mut t1[0]);
    crate::fp_generic::fpneg(&mut a[1]);

    let temp = a[0];
    fpmul_mont(&temp, &t1[0], &mut a[0]);

    let temp = a[1];
    fpmul_mont(&temp, &t1[0], &mut a[1]);
}

pub fn to_fp2mont(a: &f2elm_t, mc: &mut f2elm_t) {
    to_mont(&a[0], &mut mc[0]);
    to_mont(&a[1], &mut mc[1]);
}

pub fn from_fp2mont(ma: &f2elm_t, c: &mut f2elm_t) {
    from_mont(&ma[0], &mut c[0]);
    from_mont(&ma[1], &mut c[1]);
}

pub fn mp_add(a: &[u64], b: &[u64], c: &mut [u64], nwords: usize) -> u32 {
    let mut carry = 0;

    for i in 0..nwords {
        addc(carry, a[i], b[i], &mut carry, &mut c[i]);
    }

    carry
}

pub fn mp_shiftleft(x: &mut [u64], mut shift: u32, nwords: usize) {
    let mut j: u32 = 0;

    while shift > crate::RADIX as u32 {
        j += 1;
        shift -= crate::RADIX as u32;
    }

    for i in 0..nwords - j as usize {
        x[nwords - 1 - i] = x[nwords - 1 - i - j as usize];
    }

    for i in nwords - j as usize..nwords {
        x[nwords - 1 - i] = 0;
    }

    if shift != 0 {
        let mut j = nwords - 1;

        while j > 0 {
            shiftl(x[j], x[j - 1], shift, &mut x[j]);
            j -= 1;
        }

        x[0] <<= shift;
    }
}

pub fn mp_shiftr1(x: &mut [u64], nwords: usize) {
    for i in 0..nwords.wrapping_sub(1) {
        x[i] = x[i] >> 1 as i32 ^ x[i + 1] << crate::RADIX - 1;
    }

    x[nwords - 1] >>= 1 as i32;
}

pub fn mp_shiftl1(x: &mut [u64], nwords: usize) {
    let mut i = nwords - 1;

    while i > 0 {
        x[i] = x[i] << 1 as i32 ^ x[i - 1] >> (crate::RADIX - 1) as i32;
        i -= 1
    }

    x[0] <<= 1 as i32;
}

#[cfg(test)]
pub mod tests {
    pub type uint64_t = libc::c_ulong;
    pub type digit_t = uint64_t;
    type felm_t = [u64; 12];
    type f2elm_t = [felm_t; 2];
    type dfelm_t = [u64; 24];

    use crate::fp_generic::tests::*;
    use crate::util::fillbytes_u64;

    #[test]
    fn test_fpx() {
        for _ in 0..500 {
            unsafe {
                let mut a = [0u64; 12];
                let mut mc = [0u64; 12];
                let mut md = [0u64; 12];

                fillbytes_u64(&mut [&mut a]);

                crate::fpx::to_mont(&a, &mut mc);
                to_mont(a.as_ptr(), md.as_mut_ptr());
            }

            unsafe {
                let mut ma = [0u64; 12];
                let mut c = [0u64; 12];
                let mut d = [0u64; 12];

                fillbytes_u64(&mut [&mut ma]);

                crate::fpx::from_mont(&ma, &mut c);
                from_mont(ma.as_ptr(), d.as_mut_ptr());
            }

            unsafe {
                let mut ma = [0u64; 12];
                let mut mb = [0u64; 12];
                let mut mc = [0u64; 12];
                let mut md = [0u64; 12];

                fillbytes_u64(&mut [&mut ma, &mut mb]);

                crate::fpx::fpmul_mont(&ma, &mb, &mut mc);
                fpmul_mont(ma.as_ptr(), mb.as_ptr(), md.as_mut_ptr());
            }

            unsafe {
                let mut x = [0u64; 100];
                fillbytes_u64(&mut [&mut x]);

                let mut y = x;

                crate::fpx::mp_shiftleft(&mut x, 12, 50);
                mp_shiftleft(y.as_mut_ptr(), 12, 50);

                assert_eq!(x.as_ref(), y.as_ref());
            }

            unsafe {
                let mut x = [0u64; 100];
                fillbytes_u64(&mut [&mut x]);

                let mut y = x;

                crate::fpx::mp_shiftr1(&mut x, 50);
                mp_shiftr1(y.as_mut_ptr(), 50);

                assert_eq!(x.as_ref(), y.as_ref());
            }

            unsafe {
                let mut x = [0u64; 100];
                fillbytes_u64(&mut [&mut x]);

                let mut y = x;

                crate::fpx::mp_shiftl1(&mut x, 50);
                mp_shiftl1(y.as_mut_ptr(), 50);

                assert_eq!(x.as_ref(), y.as_ref());
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpcopy(mut a: *const digit_t, mut c: *mut digit_t) {
        // Copy a field element, c = a.
        let mut i: libc::c_uint = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            *c.offset(i as isize) = *a.offset(i as isize);
            i = i.wrapping_add(1)
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpzero(mut a: *mut digit_t) {
        // Zero a field element, a = 0.
        let mut i: libc::c_uint = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            *a.offset(i as isize) = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn copy_words(
        mut a: *const digit_t,
        mut c: *mut digit_t,
        nwords: libc::c_uint,
    ) {
        // Copy wordsize digits, c = a, where lng(a) = nwords.
        let mut i: libc::c_uint = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < nwords {
            *c.offset(i as isize) = *a.offset(i as isize);
            i = i.wrapping_add(1)
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpmul_mont(
        mut ma: *const digit_t,
        mut mb: *const digit_t,
        mut mc: *mut digit_t,
    ) {
        // Multiprecision multiplication, c = a*b mod p.
        let mut temp: dfelm_t = [
            0 as libc::c_int as digit_t,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ];
        mp_mul(ma, mb, temp.as_mut_ptr(), 12 as libc::c_int as libc::c_uint);
        rdc_mont(temp.as_mut_ptr(), mc);
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpsqr_mont(mut ma: *const digit_t, mut mc: *mut digit_t) {
        // Multiprecision squaring, c = a^2 mod p.
        let mut temp: dfelm_t = [
            0 as libc::c_int as digit_t,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ];
        mp_mul(ma, ma, temp.as_mut_ptr(), 12 as libc::c_int as libc::c_uint);
        rdc_mont(temp.as_mut_ptr(), mc);
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2copy(mut a: *const felm_t, mut c: *mut felm_t) {
        // Copy a GF(p^2) element, c = a.
        fpcopy(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpcopy(
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2zero(mut a: *mut felm_t) {
        // Zero a GF(p^2) element, a = 0.
        fpzero((*a.offset(0 as libc::c_int as isize)).as_mut_ptr());
        fpzero((*a.offset(1 as libc::c_int as isize)).as_mut_ptr());
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2neg(mut a: *mut felm_t) {
        // GF(p^2) negation, a = -a in GF(p^2).
        fpneg((*a.offset(0 as libc::c_int as isize)).as_mut_ptr());
        fpneg((*a.offset(1 as libc::c_int as isize)).as_mut_ptr());
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2add(
        mut a: *const felm_t,
        mut b: *const felm_t,
        mut c: *mut felm_t,
    ) {
        // GF(p^2) addition, c = a+b in GF(p^2).
        fpadd(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*b.offset(0 as libc::c_int as isize)).as_ptr(),
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpadd(
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*b.offset(1 as libc::c_int as isize)).as_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2sub(
        mut a: *const felm_t,
        mut b: *const felm_t,
        mut c: *mut felm_t,
    ) {
        // GF(p^2) subtraction, c = a-b in GF(p^2).
        fpsub(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*b.offset(0 as libc::c_int as isize)).as_ptr(),
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpsub(
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*b.offset(1 as libc::c_int as isize)).as_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2div2(mut a: *const felm_t, mut c: *mut felm_t) {
        // GF(p^2) division by two, c = a/2  in GF(p^2).

        fpdiv2(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpdiv2(
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2correction(mut a: *mut felm_t) {
        // Modular correction, a = a in GF(p^2).
        fpcorrection((*a.offset(0 as libc::c_int as isize)).as_mut_ptr());
        fpcorrection((*a.offset(1 as libc::c_int as isize)).as_mut_ptr());
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_add(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
        nwords: libc::c_uint,
    ) -> libc::c_uint {
        // Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit.
        let mut i: libc::c_uint = 0;
        let mut carry: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < nwords {
            let mut tempReg: uint128_t = (*a.offset(i as isize) as uint128_t)
                .wrapping_add(*b.offset(i as isize) as uint128_t)
                .wrapping_add(carry as uint128_t);
            carry = (tempReg >> 64 as libc::c_int) as digit_t as libc::c_uint;
            *c.offset(i as isize) = tempReg as digit_t;
            i = i.wrapping_add(1)
        }
        return carry;
    }

    unsafe extern "C" fn mp_addfast(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
    ) {
        // Multiprecision addition, c = a+b.
        mp_add(a, b, c, 12 as libc::c_int as libc::c_uint);
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2sqr_mont(mut a: *const felm_t, mut c: *mut felm_t) {
        // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
        // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1]
        // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1]
        let mut t1: felm_t = [0; 12]; // t1 = a0+a1
        let mut t2: felm_t = [0; 12]; // t2 = a0-a1
        let mut t3: felm_t = [0; 12]; // t3 = 2a0
        mp_addfast(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            t1.as_mut_ptr(),
        ); // c0 = (a0+a1)(a0-a1)
        fpsub(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            t2.as_mut_ptr(),
        );
        mp_addfast(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            t3.as_mut_ptr(),
        );
        fpmul_mont(
            t1.as_mut_ptr() as *const digit_t,
            t2.as_mut_ptr() as *const digit_t,
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpmul_mont(
            t3.as_mut_ptr() as *const digit_t,
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
        // c1 = 2a0*a1
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_sub(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
        nwords: libc::c_uint,
    ) -> libc::c_uint {
        // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit.
        let mut i: libc::c_uint = 0;
        let mut borrow: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < nwords {
            let mut tempReg: uint128_t = (*a.offset(i as isize) as uint128_t)
                .wrapping_sub(*b.offset(i as isize) as uint128_t)
                .wrapping_sub(borrow as uint128_t);
            borrow = (tempReg
                >> (::std::mem::size_of::<uint128_t>() as libc::c_ulong)
                    .wrapping_mul(8 as libc::c_int as libc::c_ulong)
                    .wrapping_sub(1 as libc::c_int as libc::c_ulong))
                as digit_t as libc::c_uint;
            *c.offset(i as isize) = tempReg as digit_t;
            i = i.wrapping_add(1)
        }
        return borrow;
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_subfast(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
    ) -> digit_t {
        // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD.
        // If c < 0 then returns mask = 0xFF..F, else mask = 0x00..0
        return (0 as libc::c_int as libc::c_ulong).wrapping_sub(mp_sub(
            a,
            b,
            c,
            (2 as libc::c_int * 12 as libc::c_int) as libc::c_uint,
        ) as digit_t);
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_dblsubfast(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
    ) {
        // Multiprecision subtraction, c = c-a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD.
        // Inputs should be s.t. c > a and c > b
        mp_sub(
            c,
            a,
            c,
            (2 as libc::c_int * 12 as libc::c_int) as libc::c_uint,
        );
        mp_sub(
            c,
            b,
            c,
            (2 as libc::c_int * 12 as libc::c_int) as libc::c_uint,
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2mul_mont(
        mut a: *const felm_t,
        mut b: *const felm_t,
        mut c: *mut felm_t,
    ) {
        // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
        // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1]
        // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1]
        let mut t1: felm_t = [0; 12]; // t1 = a0+a1
        let mut t2: felm_t = [0; 12]; // t2 = b0+b1
        let mut tt1: dfelm_t = [0; 24]; // tt1 = a0*b0
        let mut tt2: dfelm_t = [0; 24]; // tt2 = a1*b1
        let mut tt3: dfelm_t = [0; 24]; // tt3 = (a0+a1)*(b0+b1)
        let mut mask: digit_t = 0; // tt3 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1
        let mut i: libc::c_uint = 0; // tt1 = a0*b0 - a1*b1. If tt1 < 0 then mask = 0xFF..F, else if tt1 >= 0 then mask = 0x00..0
        mp_addfast(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            t1.as_mut_ptr(),
        ); // c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1
        mp_addfast(
            (*b.offset(0 as libc::c_int as isize)).as_ptr(),
            (*b.offset(1 as libc::c_int as isize)).as_ptr(),
            t2.as_mut_ptr(),
        );
        mp_mul(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*b.offset(0 as libc::c_int as isize)).as_ptr(),
            tt1.as_mut_ptr(),
            12 as libc::c_int as libc::c_uint,
        );
        mp_mul(
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*b.offset(1 as libc::c_int as isize)).as_ptr(),
            tt2.as_mut_ptr(),
            12 as libc::c_int as libc::c_uint,
        );
        mp_mul(
            t1.as_mut_ptr(),
            t2.as_mut_ptr(),
            tt3.as_mut_ptr(),
            12 as libc::c_int as libc::c_uint,
        );
        mp_dblsubfast(tt1.as_mut_ptr(), tt2.as_mut_ptr(), tt3.as_mut_ptr());
        mask = mp_subfast(tt1.as_mut_ptr(), tt2.as_mut_ptr(), tt1.as_mut_ptr());
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            t1[i as usize] = *(crate::p751.as_ptr() as *mut digit_t).offset(i as isize) & mask;
            i = i.wrapping_add(1)
        }
        rdc_mont(
            tt3.as_mut_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
        mp_addfast(
            &mut *tt1.as_mut_ptr().offset(12 as libc::c_int as isize) as *mut digit_t,
            t1.as_mut_ptr(),
            &mut *tt1.as_mut_ptr().offset(12 as libc::c_int as isize) as *mut digit_t,
        );
        rdc_mont(
            tt1.as_mut_ptr(),
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        // c[0] = a0*b0 - a1*b1
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpinv_chain_mont(mut a: *mut digit_t) {
        // Chain to compute a^(p-3)/4 using Montgomery arithmetic.
        let mut i: libc::c_uint = 0;
        let mut j: libc::c_uint = 0;
        let mut t: [felm_t; 27] = [[0; 12]; 27];
        let mut tt: felm_t = [0; 12];
        // Precomputed table
        fpsqr_mont(a as *const digit_t, tt.as_mut_ptr());
        fpmul_mont(
            a as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[0 as libc::c_int as usize].as_mut_ptr(),
        );
        fpmul_mont(
            t[0 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[1 as libc::c_int as usize].as_mut_ptr(),
        );
        fpmul_mont(
            t[1 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[2 as libc::c_int as usize].as_mut_ptr(),
        );
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[3 as libc::c_int as usize].as_mut_ptr(),
        );
        fpmul_mont(
            t[3 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[3 as libc::c_int as usize].as_mut_ptr(),
        );
        i = 3 as libc::c_int as libc::c_uint;
        while i <= 8 as libc::c_int as libc::c_uint {
            //
            fpmul_mont(
                t[i as usize].as_mut_ptr() as *const digit_t,
                tt.as_mut_ptr() as *const digit_t,
                t[i.wrapping_add(1 as libc::c_int as libc::c_uint) as usize].as_mut_ptr(),
            );
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[9 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[9 as libc::c_int as usize].as_mut_ptr(),
        );
        i = 9 as libc::c_int as libc::c_uint;
        while i <= 20 as libc::c_int as libc::c_uint {
            fpmul_mont(
                t[i as usize].as_mut_ptr() as *const digit_t,
                tt.as_mut_ptr() as *const digit_t,
                t[i.wrapping_add(1 as libc::c_int as libc::c_uint) as usize].as_mut_ptr(),
            );
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[21 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[21 as libc::c_int as usize].as_mut_ptr(),
        );
        i = 21 as libc::c_int as libc::c_uint;
        while i <= 24 as libc::c_int as libc::c_uint {
            fpmul_mont(
                t[i as usize].as_mut_ptr() as *const digit_t,
                tt.as_mut_ptr() as *const digit_t,
                t[i.wrapping_add(1 as libc::c_int as libc::c_uint) as usize].as_mut_ptr(),
            );
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[25 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[25 as libc::c_int as usize].as_mut_ptr(),
        );
        fpmul_mont(
            t[25 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            t[26 as libc::c_int as usize].as_mut_ptr(),
        );
        fpcopy(a as *const digit_t, tt.as_mut_ptr());
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[20 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[24 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[11 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[8 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[23 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 9 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 10 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[15 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[13 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[26 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[20 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[11 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[10 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[14 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[4 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 10 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[18 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[1 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[22 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 10 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[6 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[24 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[9 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[18 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[17 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            a as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 10 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[16 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[7 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[0 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[12 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[19 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[22 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[25 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[10 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[22 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[18 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[4 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[14 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[13 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[5 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[23 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[21 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[23 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[12 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[9 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[3 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[13 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[17 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[26 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[5 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 8 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[8 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[2 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 6 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[11 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        i = 0 as libc::c_int as libc::c_uint;
        while i < 7 as libc::c_int as libc::c_uint {
            fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
            i = i.wrapping_add(1)
        }
        fpmul_mont(
            t[20 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr() as *const digit_t,
            tt.as_mut_ptr(),
        );
        j = 0 as libc::c_int as libc::c_uint;
        while j < 61 as libc::c_int as libc::c_uint {
            i = 0 as libc::c_int as libc::c_uint;
            while i < 6 as libc::c_int as libc::c_uint {
                fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
                i = i.wrapping_add(1)
            }
            fpmul_mont(
                t[26 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
                tt.as_mut_ptr() as *const digit_t,
                tt.as_mut_ptr(),
            );
            j = j.wrapping_add(1)
        }
        fpcopy(tt.as_mut_ptr() as *const digit_t, a);
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpinv_mont(mut a: *mut digit_t) {
        // Field inversion using Montgomery arithmetic, a = a^(-1)*R mod p.
        let mut tt: felm_t = [0; 12];
        fpcopy(a as *const digit_t, tt.as_mut_ptr());
        fpinv_chain_mont(tt.as_mut_ptr());
        fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
        fpsqr_mont(tt.as_mut_ptr() as *const digit_t, tt.as_mut_ptr());
        fpmul_mont(a as *const digit_t, tt.as_mut_ptr() as *const digit_t, a);
    }

    #[no_mangle]
    pub unsafe extern "C" fn fp2inv_mont(mut a: *mut felm_t) {
        // GF(p^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2).
        let mut t1: f2elm_t = [[0; 12]; 2]; // t10 = a0^2
        fpsqr_mont(
            (*a.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const digit_t,
            t1[0 as libc::c_int as usize].as_mut_ptr(),
        ); // t11 = a1^2
        fpsqr_mont(
            (*a.offset(1 as libc::c_int as isize)).as_mut_ptr() as *const digit_t,
            t1[1 as libc::c_int as usize].as_mut_ptr(),
        ); // t10 = a0^2+a1^2
        fpadd(
            t1[0 as libc::c_int as usize].as_mut_ptr(),
            t1[1 as libc::c_int as usize].as_mut_ptr(),
            t1[0 as libc::c_int as usize].as_mut_ptr(),
        ); // t10 = (a0^2+a1^2)^-1
        fpinv_mont(t1[0 as libc::c_int as usize].as_mut_ptr()); // a = a0-i*a1
        fpneg((*a.offset(1 as libc::c_int as isize)).as_mut_ptr());
        fpmul_mont(
            (*a.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const digit_t,
            t1[0 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            (*a.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        fpmul_mont(
            (*a.offset(1 as libc::c_int as isize)).as_mut_ptr() as *const digit_t,
            t1[0 as libc::c_int as usize].as_mut_ptr() as *const digit_t,
            (*a.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
        // a = (a0-i*a1)*(a0^2+a1^2)^-1
    }

    #[no_mangle]
    pub unsafe extern "C" fn to_mont(mut a: *const digit_t, mut mc: *mut digit_t) {
        // Conversion to Montgomery representation,
        // mc = a*R^2*R^(-1) mod p = a*R mod p, where a in [0, p-1].
        // The Montgomery constant R^2 mod p is the global value "Montgomery_R2".
        fpmul_mont(
            a,
            &crate::Montgomery_R2 as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            mc,
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn from_mont(mut ma: *const digit_t, mut c: *mut digit_t) {
        // Conversion from Montgomery representation to standard representation,
        // c = ma*R^(-1) mod p = a mod p, where ma in [0, p-1].
        let mut one: [digit_t; 12] = [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        one[0 as libc::c_int as usize] = 1 as libc::c_int as digit_t;
        fpmul_mont(ma, one.as_mut_ptr() as *const digit_t, c);
        fpcorrection(c);
    }

    #[no_mangle]
    pub unsafe extern "C" fn to_fp2mont(mut a: *const felm_t, mut mc: *mut felm_t) {
        // Conversion of a GF(p^2) element to Montgomery representation,
        // mc_i = a_i*R^2*R^(-1) = a_i*R in GF(p^2).
        to_mont(
            (*a.offset(0 as libc::c_int as isize)).as_ptr(),
            (*mc.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        to_mont(
            (*a.offset(1 as libc::c_int as isize)).as_ptr(),
            (*mc.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn from_fp2mont(mut ma: *const felm_t, mut c: *mut felm_t) {
        // Conversion of a GF(p^2) element from Montgomery representation to standard representation,
        // c_i = ma_i*R^(-1) = a_i in GF(p^2).
        from_mont(
            (*ma.offset(0 as libc::c_int as isize)).as_ptr(),
            (*c.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        );
        from_mont(
            (*ma.offset(1 as libc::c_int as isize)).as_ptr(),
            (*c.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        );
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_shiftleft(
        mut x: *mut digit_t,
        mut shift: libc::c_uint,
        nwords: libc::c_uint,
    ) {
        let mut i: libc::c_uint = 0;
        let mut j: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        while shift > 64 as libc::c_int as libc::c_uint {
            j = j.wrapping_add(1 as libc::c_int as libc::c_uint);
            shift = shift.wrapping_sub(64 as libc::c_int as libc::c_uint)
        }
        i = 0 as libc::c_int as libc::c_uint;
        while i < nwords.wrapping_sub(j) {
            *x.offset(
                nwords
                    .wrapping_sub(1 as libc::c_int as libc::c_uint)
                    .wrapping_sub(i) as isize,
            ) = *x.offset(
                nwords
                    .wrapping_sub(1 as libc::c_int as libc::c_uint)
                    .wrapping_sub(i)
                    .wrapping_sub(j) as isize,
            );
            i = i.wrapping_add(1)
        }
        i = nwords.wrapping_sub(j);
        while i < nwords {
            *x.offset(
                nwords
                    .wrapping_sub(1 as libc::c_int as libc::c_uint)
                    .wrapping_sub(i) as isize,
            ) = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
        if shift != 0 as libc::c_int as libc::c_uint {
            j = nwords.wrapping_sub(1 as libc::c_int as libc::c_uint);
            while j > 0 as libc::c_int as libc::c_uint {
                *x.offset(j as isize) = *x.offset(j as isize) << shift
                    ^ *x.offset(j.wrapping_sub(1 as libc::c_int as libc::c_uint) as isize)
                        >> (64 as libc::c_int as libc::c_uint).wrapping_sub(shift);
                j = j.wrapping_sub(1)
            }
            *x.offset(0 as libc::c_int as isize) <<= shift
        };
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_shiftr1(mut x: *mut digit_t, nwords: libc::c_uint) {
        // Multiprecision right shift by one.
        let mut i: libc::c_uint = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < nwords.wrapping_sub(1 as libc::c_int as libc::c_uint) {
            *x.offset(i as isize) = *x.offset(i as isize) >> 1 as libc::c_int
                ^ *x.offset(i.wrapping_add(1 as libc::c_int as libc::c_uint) as isize)
                    << 64 as libc::c_int - 1 as libc::c_int;
            i = i.wrapping_add(1)
        }
        *x.offset(nwords.wrapping_sub(1 as libc::c_int as libc::c_uint) as isize) >>=
            1 as libc::c_int;
    }

    #[no_mangle]
    pub unsafe extern "C" fn mp_shiftl1(mut x: *mut digit_t, nwords: libc::c_uint) {
        // Multiprecision left shift by one.
        let mut i: libc::c_int = 0;
        i = nwords.wrapping_sub(1 as libc::c_int as libc::c_uint) as libc::c_int;
        while i > 0 as libc::c_int {
            *x.offset(i as isize) = *x.offset(i as isize) << 1 as libc::c_int
                ^ *x.offset((i - 1 as libc::c_int) as isize)
                    >> 64 as libc::c_int - 1 as libc::c_int;
            i -= 1
        }
        *x.offset(0 as libc::c_int as isize) <<= 1 as libc::c_int;
    }
}
