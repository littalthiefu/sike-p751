use crate::util::{addc, shiftl, subc};
use crate::{f2elm_t, felm_t};

pub fn to_mont(a: &felm_t, mc: &mut felm_t) {
    fpmul_mont(a, &crate::Montgomery_R2, mc);
}

pub fn from_mont(ma: &felm_t, c: &mut felm_t) {
    let mut one = [0u64; crate::NWORDS_FIELD];
    one[0] = 1;

    fpmul_mont(ma, &one, c);
    crate::fp_generic::fpcorrection751(c);
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
    crate::fp_generic::fpneg751(&mut a[0]);
    crate::fp_generic::fpneg751(&mut a[1]);
}

pub fn fp2add(a: &f2elm_t, b: &f2elm_t, c: &mut f2elm_t) {
    crate::fp_generic::fpadd751(&a[0], &b[0], &mut c[0]);
    crate::fp_generic::fpadd751(&a[1], &b[1], &mut c[1]);
}

pub fn fp2sub(a: &f2elm_t, b: &f2elm_t, c: &mut f2elm_t) {
    crate::fp_generic::fpsub751(&a[0], &b[0], &mut c[0]);
    crate::fp_generic::fpsub751(&a[1], &b[1], &mut c[1]);
}

pub fn fp2div2(a: &f2elm_t, c: &mut f2elm_t) {
    crate::fp_generic::fpdiv2_751(&a[0], &mut c[0]);
}

pub fn fp2correction(a: &mut f2elm_t) {
    crate::fp_generic::fpcorrection751(&mut a[0]);
    crate::fp_generic::fpcorrection751(&mut a[1]);
}

pub fn mp_addfast(a: &[u64], b: &[u64], c: &mut [u64]) {
    mp_add(a, b, c, 12);
}

pub fn fp2sqr_mont(a: &f2elm_t, c: &mut f2elm_t) {
    let mut t1 = [0u64; 12];
    let mut t2 = [0u64; 12];
    let mut t3 = [0u64; 12];

    mp_addfast(&a[0], &a[1], &mut t1);
    crate::fp_generic::fpsub751(&a[0], &a[1], &mut t2);
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
    crate::fp_generic::fpadd751(&temp[0], &temp[1], &mut t1[0]);
    fpinv_mont(&mut t1[0]);
    crate::fp_generic::fpneg751(&mut a[1]);

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
mod tests {
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
