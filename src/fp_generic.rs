use crate::util::{addc, mul, subc};

pub fn fpaddx(a: *const u64, b: *const u64, c: &mut [u64]) {
    unsafe {
        let mut carry: u32 = 0;

        for i in 0..12 {
            addc(
                carry,
                *a.offset(i),
                *b.offset(i),
                &mut carry,
                &mut c[i as usize],
            );
        }

        carry = 0;

        for i in 0..12 {
            subc(
                carry,
                c[i],
                crate::p751x2[i as usize],
                &mut carry,
                &mut c[i],
            );
        }

        let mask = (0 as u64).wrapping_sub(carry as u64);

        carry = 0;

        for i in 0..12 {
            addc(
                carry,
                c[i],
                crate::p751x2[i as usize] & mask,
                &mut carry,
                &mut c[i],
            );
        }
    }
}

pub fn fpadd(a: &[u64], b: &[u64], c: &mut [u64]) {
    let mut carry: u32 = 0;

    for i in 0..12 {
        addc(carry, a[i], b[i], &mut carry, &mut c[i]);
    }

    carry = 0;

    for i in 0..12 {
        subc(carry, c[i], crate::p751x2[i], &mut carry, &mut c[i]);
    }

    let mask = (0 as u64).wrapping_sub(carry as u64);

    carry = 0;

    for i in 0..12 {
        addc(carry, c[i], crate::p751x2[i] & mask, &mut carry, &mut c[i]);
    }
}

pub fn fpsubx(a: *const u64, b: *const u64, c: &mut [u64]) {
    unsafe {
        let mut borrow: u32 = 0;

        for i in 0..12 {
            subc(
                borrow,
                *a.offset(i),
                *b.offset(i),
                &mut borrow,
                &mut c[i as usize],
            );
        }

        let mask = (0 as u64).wrapping_sub(borrow as u64);
        borrow = 0;

        for i in 0..12 {
            addc(
                borrow,
                c[i],
                crate::p751x2[i] & mask,
                &mut borrow,
                &mut c[i],
            );
        }
    }
}

pub fn fpsub(a: &[u64], b: &[u64], c: &mut [u64]) {
    let mut borrow: u32 = 0;

    for i in 0..12 {
        subc(borrow, a[i], b[i], &mut borrow, &mut c[i]);
    }

    let mask = (0 as u64).wrapping_sub(borrow as u64);
    borrow = 0;

    for i in 0..12 {
        addc(
            borrow,
            c[i],
            crate::p751x2[i] & mask,
            &mut borrow,
            &mut c[i],
        );
    }
}

pub fn fpneg(a: &mut [u64]) {
    let mut borrow: u32 = 0;

    for i in 0..12 {
        subc(borrow, crate::p751x2[i], a[i], &mut borrow, &mut a[i]);
    }
}

pub fn fpdiv2(a: &[u64], c: &mut [u64]) {
    let mut carry: u32 = 0;
    let mask = (0 as u64).wrapping_sub(a[0] & 1 as u64);

    for i in 0..crate::NWORDS_FIELD {
        addc(carry, a[i], crate::p751[i] & mask, &mut carry, &mut c[i]);
    }

    crate::fpx::mp_shiftr1(c, 12);
}

pub fn fpcorrection(a: &mut [u64]) {
    let mut borrow: u32 = 0;

    for i in 0..12 {
        let temp_reg = (a[i] as u128)
            .wrapping_sub(crate::p751[i] as u128)
            .wrapping_sub(borrow as u128);

        borrow = (temp_reg >> (16 as u64).wrapping_mul(8).wrapping_sub(1)) as u32;
        a[i] = temp_reg as u64;
    }

    let mask = (0 as u64).wrapping_sub(borrow as u64);
    borrow = 0;

    for i in 0..12 {
        let temp_reg = (a[i] as u128)
            .wrapping_add((crate::p751[i] & mask) as u128)
            .wrapping_add(borrow as u128);

        borrow = (temp_reg >> 64 as i32) as u32;
        a[i] = temp_reg as u64;
    }
}

pub fn digit_x_digit(a: u64, b: u64, c: &mut [u64]) {
    let mut al: u64 = 0;
    let mut ah: u64 = 0;
    let mut bl: u64 = 0;
    let mut bh: u64 = 0;
    let mut temp: u64 = 0;
    let mut albl: u64 = 0;
    let mut albh: u64 = 0;
    let mut ahbl: u64 = 0;
    let mut ahbh: u64 = 0;
    let mut res1: u64 = 0;
    let mut res2: u64 = 0;
    let mut res3: u64 = 0;
    let mut carry: u64 = 0;
    let mut mask_low: u64 = -(1 as i32) as u64 >> (8 as u64).wrapping_mul(4);
    let mut mask_high: u64 = (-(1 as i32) as u64) << (8 as u64).wrapping_mul(4);

    al = a & mask_low;
    ah = a >> (8 as u64).wrapping_mul(4);
    bl = b & mask_low;
    bh = b >> (8 as u64).wrapping_mul(4);

    albl = al.wrapping_mul(bl);
    albh = al.wrapping_mul(bh);
    ahbl = ah.wrapping_mul(bl);
    ahbh = ah.wrapping_mul(bh);
    c[0] = albl & mask_low;

    res1 = albl >> (8 as u64).wrapping_mul(4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;
    temp = res1.wrapping_add(res2).wrapping_add(res3);
    carry = temp >> (8 as u64).wrapping_mul(4);
    c[0] ^= temp << (8 as u64).wrapping_mul(4);

    res1 = ahbl >> (8 as u64).wrapping_mul(4);
    res2 = albh >> (8 as u64).wrapping_mul(4);
    res3 = ahbh & mask_low;
    temp = res1
        .wrapping_add(res2)
        .wrapping_add(res3)
        .wrapping_add(carry);

    c[1] = temp & mask_low;
    carry = temp & mask_high;
    c[1] ^= (ahbh & mask_high).wrapping_add(carry);
}

pub fn mp_mul(a: &[u64], b: &[u64], c: &mut [u64], nwords: u32) {
    let mut t: u64 = 0;
    let mut u: u64 = 0;
    let mut v: u64 = 0;
    let mut uv: [u64; 2] = [0; 2];
    let mut carry: u32 = 0;

    for i in 0..nwords {
        for j in 0..i + 1 {
            let mut temp = 0;

            mul(
                a[j as usize],
                b[i.wrapping_sub(j) as usize],
                &mut uv[1],
                &mut temp,
            );

            uv[0] = temp;

            addc(0, uv[0], v, &mut carry, &mut v);
            addc(carry, uv[1], u, &mut carry, &mut u);

            t = t.wrapping_add(carry as u64);
        }

        c[i as usize] = v;
        v = u;
        u = t;
        t = 0;
    }

    for i in nwords..(2 as u32).wrapping_mul(nwords).wrapping_sub(1) {
        for j in i.wrapping_sub(nwords).wrapping_add(1)..nwords {
            let mut temp = 0;

            mul(
                a[j as usize],
                b[i.wrapping_sub(j) as usize],
                &mut uv[1],
                &mut temp,
            );

            uv[0] = temp;

            addc(0, uv[0], v, &mut carry, &mut v);
            addc(carry, uv[1], u, &mut carry, &mut u);

            t = t.wrapping_add(carry as u64);
        }

        c[i as usize] = v;
        v = u;
        u = t;
        t = 0;
    }

    c[(2 as u32).wrapping_mul(nwords).wrapping_sub(1) as usize] = v;
}

pub fn rdc_mont(ma: &[u64], mc: &mut [u64]) {
    let mut carry: u32 = 0;
    let mut count: u32 = 5;
    let mut UV: [u64; 2] = [0; 2];
    let mut t: u64 = 0;
    let mut u: u64 = 0;
    let mut v: u64 = 0;

    for i in 0..12 {
        mc[i as usize] = 0;
    }

    for i in 0..12 {
        for j in 0..i {
            if j < (i as u32).wrapping_sub(5).wrapping_add(1) {
                let mut temp = 0;
                mul(
                    mc[j as usize],
                    crate::p751p1[i.wrapping_sub(j) as usize],
                    &mut UV[1],
                    &mut temp,
                );

                UV[0] = temp;

                addc(0, UV[0], v, &mut carry, &mut v);
                addc(carry, UV[1], u, &mut carry, &mut u);
                t = (t as u64).wrapping_add(carry as u64);
            }
        }

        addc(0, v, ma[i as usize], &mut carry, &mut v);
        addc(carry, u, 0, &mut carry, &mut u);

        t = (t as u64).wrapping_add(carry as u64);
        mc[i as usize] = v;
        v = u;
        u = t;
        t = 0;
    }

    for i in 12..23 as u32 {
        if count > 0 {
            count = count.wrapping_sub(1 as i32 as u32)
        }

        for j in i.wrapping_sub(12).wrapping_add(1)..12 {
            if j < (12 as u32).wrapping_sub(count) {
                let mut temp = 0;
                mul(
                    mc[j as usize],
                    crate::p751p1[i.wrapping_sub(j) as usize],
                    &mut UV[1],
                    &mut temp,
                );

                UV[0] = temp;

                addc(0, UV[0], v, &mut carry, &mut v);
                addc(carry, UV[1], u, &mut carry, &mut u);
                t = (t as u64).wrapping_add(carry as u64);
            }
        }

        addc(0, v, ma[i as usize], &mut carry, &mut v);
        addc(carry, u, 0, &mut carry, &mut u);

        t = (t as u64).wrapping_add(carry as u64);
        mc[(i - 12) as usize] = v;
        v = u;
        u = t;
        t = 0;
    }

    addc(0, v, ma[23], &mut carry, &mut v);
    mc[11] = v;
}

#[cfg(test)]
pub mod tests {
    use crate::*;

    #[test]
    fn test_fp_generic() {
        for _ in 0..500 {
            let mut a = [0u64; 12];
            let mut b = [0u64; 12];

            let mut c = [0u64; 12];
            let mut d = [0u64; 12];

            crate::util::fillbytes_u64(&mut [&mut a, &mut b]);

            unsafe {
                crate::fp_generic::fpadd(&a, &b, &mut c);
                fpadd(a.as_ptr(), b.as_ptr(), d.as_mut_ptr());
                assert_eq!(c, d);
            }

            unsafe {
                crate::fp_generic::fpsub(&a, &b, &mut c);
                fpsub(a.as_ptr(), b.as_ptr(), d.as_mut_ptr());
                assert_eq!(c, d);
            }

            unsafe {
                c = a;
                d = a;

                crate::fp_generic::fpneg(&mut c);
                fpneg(d.as_mut_ptr());
                assert_eq!(c, d);
            }

            unsafe {
                let mut a = [0u64; 12];
                let mut b = [0u64; 12];
                let mut c = [0u64; 12];

                crate::util::fillbytes_u64(&mut [&mut a]);

                crate::fp_generic::fpdiv2(&a, &mut b);
                fpdiv2(a.as_ptr(), c.as_mut_ptr());

                assert_eq!(b, c);
            }

            unsafe {
                c = a;
                d = a;

                crate::fp_generic::fpcorrection(&mut c);
                fpcorrection(d.as_mut_ptr());
                assert_eq!(c, d);
            }

            unsafe {
                let mut c = [0u64; 2];
                let mut d = [0u64; 2];

                crate::fp_generic::digit_x_digit(a[0], b[0], &mut c);
                digit_x_digit(a[0], b[0], d.as_mut_ptr());

                assert_eq!(c, d);
            }

            unsafe {
                let mut a = [0u64; 20];
                let mut b = [0u64; 20];

                let mut c = [0u64; 40];
                let mut d = [0u64; 40];

                crate::util::fillbytes_u64(&mut [&mut a, &mut b]);
                crate::fp_generic::mp_mul(&a, &b, &mut c, 10);
                mp_mul(a.as_ptr(), b.as_ptr(), d.as_mut_ptr(), 10);

                assert_eq!(c.as_ref(), d.as_ref());
            }

            unsafe {
                let mut ma = [0u64; 50];
                let mut mc = [0u64; 50];
                let mut md = [0u64; 50];

                crate::util::fillbytes_u64(&mut [&mut ma]);
                crate::fp_generic::rdc_mont(&ma, &mut mc);
                rdc_mont(ma.as_ptr(), md.as_mut_ptr());

                assert_eq!(mc.as_ref(), md.as_ref());
            }
        }
    }

    pub unsafe extern "C" fn fpadd(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
    ) {
        // Modular addition, c = a+b mod p751.
        // Inputs: a, b in [0, 2*p751-1]
        // Output: c in [0, 2*p751-1]
        let mut i: libc::c_uint = 0;
        let mut carry: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        let mut mask: digit_t = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg: uint128_t = (*a.offset(i as isize) as uint128_t)
                .wrapping_add(*b.offset(i as isize) as uint128_t)
                .wrapping_add(carry as uint128_t);
            carry = (tempReg >> 64 as libc::c_int) as digit_t as libc::c_uint;
            *c.offset(i as isize) = tempReg as digit_t;
            i = i.wrapping_add(1)
        }
        carry = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg_0: uint128_t = (*c.offset(i as isize) as uint128_t)
                .wrapping_sub(
                    *(crate::p751x2.as_ptr() as *mut digit_t).offset(i as isize) as uint128_t
                )
                .wrapping_sub(carry as uint128_t);
            carry = (tempReg_0
                >> (::std::mem::size_of::<uint128_t>() as libc::c_ulong)
                    .wrapping_mul(8 as libc::c_int as libc::c_ulong)
                    .wrapping_sub(1 as libc::c_int as libc::c_ulong)) as digit_t
                as libc::c_uint;
            *c.offset(i as isize) = tempReg_0 as digit_t;
            i = i.wrapping_add(1)
        }
        mask = (0 as libc::c_int as libc::c_ulong).wrapping_sub(carry as digit_t);
        carry = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg_1: uint128_t = (*c.offset(i as isize) as uint128_t)
                .wrapping_add(
                    (*(crate::p751x2.as_ptr() as *mut digit_t).offset(i as isize) & mask)
                        as uint128_t,
                )
                .wrapping_add(carry as uint128_t);
            carry = (tempReg_1 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            *c.offset(i as isize) = tempReg_1 as digit_t;
            i = i.wrapping_add(1)
        }
    }

    pub unsafe extern "C" fn fpsub(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
    ) {
        // Modular subtraction, c = a-b mod p751.
        // Inputs: a, b in [0, 2*p751-1]
        // Output: c in [0, 2*p751-1]
        let mut i: libc::c_uint = 0;
        let mut borrow: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        let mut mask: digit_t = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
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
        mask = (0 as libc::c_int as libc::c_ulong).wrapping_sub(borrow as digit_t);
        borrow = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg_0: uint128_t = (*c.offset(i as isize) as uint128_t)
                .wrapping_add(
                    (*(crate::p751x2.as_ptr() as *mut digit_t).offset(i as isize) & mask)
                        as uint128_t,
                )
                .wrapping_add(borrow as uint128_t);
            borrow = (tempReg_0 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            *c.offset(i as isize) = tempReg_0 as digit_t;
            i = i.wrapping_add(1)
        }
    }

    pub unsafe extern "C" fn fpneg(mut a: *mut digit_t) {
        // Modular negation, a = -a mod p751.
        // Input/output: a in [0, 2*p751-1]
        let mut i: libc::c_uint = 0;
        let mut borrow: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg: uint128_t = (*(crate::p751x2.as_ptr() as *mut digit_t)
                .offset(i as isize) as uint128_t)
                .wrapping_sub(*a.offset(i as isize) as uint128_t)
                .wrapping_sub(borrow as uint128_t);
            borrow = (tempReg
                >> (::std::mem::size_of::<uint128_t>() as libc::c_ulong)
                    .wrapping_mul(8 as libc::c_int as libc::c_ulong)
                    .wrapping_sub(1 as libc::c_int as libc::c_ulong))
                as digit_t as libc::c_uint;
            *a.offset(i as isize) = tempReg as digit_t;
            i = i.wrapping_add(1)
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn fpdiv2(mut a: *const digit_t, mut c: *mut digit_t) {
        // Modular division by two, c = a/2 mod p751.
        // Input : a in [0, 2*p751-1]
        // Output: c in [0, 2*p751-1]
        let mut i: libc::c_uint = 0; // If a is odd compute a+p751
        let mut carry: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        let mut mask: digit_t = 0;
        mask = (0 as libc::c_int as libc::c_ulong)
            .wrapping_sub(*a.offset(0 as libc::c_int as isize) & 1 as libc::c_int as libc::c_ulong);
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg: uint128_t = (*a.offset(i as isize) as uint128_t)
                .wrapping_add(
                    (*(crate::p751.as_ptr() as *mut digit_t).offset(i as isize) & mask)
                        as uint128_t,
                )
                .wrapping_add(carry as uint128_t);
            carry = (tempReg >> 64 as libc::c_int) as digit_t as libc::c_uint;
            *c.offset(i as isize) = tempReg as digit_t;
            i = i.wrapping_add(1)
        }

        crate::fpx::tests::mp_shiftr1(c, 12);
    }

    pub unsafe extern "C" fn fpcorrection(mut a: *mut digit_t) {
        // Modular correction to reduce field element a in [0, 2*p751-1] to [0, p751-1].
        let mut i: libc::c_uint = 0;
        let mut borrow: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        let mut mask: digit_t = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg: uint128_t = (*a.offset(i as isize) as uint128_t)
                .wrapping_sub(
                    *(crate::p751.as_ptr() as *mut digit_t).offset(i as isize) as uint128_t
                )
                .wrapping_sub(borrow as uint128_t);
            borrow = (tempReg
                >> (::std::mem::size_of::<uint128_t>() as libc::c_ulong)
                    .wrapping_mul(8 as libc::c_int as libc::c_ulong)
                    .wrapping_sub(1 as libc::c_int as libc::c_ulong))
                as digit_t as libc::c_uint;
            *a.offset(i as isize) = tempReg as digit_t;
            i = i.wrapping_add(1)
        }
        mask = (0 as libc::c_int as libc::c_ulong).wrapping_sub(borrow as digit_t);
        borrow = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            let mut tempReg_0: uint128_t = (*a.offset(i as isize) as uint128_t)
                .wrapping_add(
                    (*(crate::p751.as_ptr() as *mut digit_t).offset(i as isize) & mask)
                        as uint128_t,
                )
                .wrapping_add(borrow as uint128_t);
            borrow = (tempReg_0 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            *a.offset(i as isize) = tempReg_0 as digit_t;
            i = i.wrapping_add(1)
        }
    }

    pub unsafe extern "C" fn digit_x_digit(a: digit_t, b: digit_t, mut c: *mut digit_t) {
        // Digit multiplication, digit * digit -> 2-digit result
        let mut al: digit_t = 0; // Low part
        let mut ah: digit_t = 0; // High part
        let mut bl: digit_t = 0; // C00
        let mut bh: digit_t = 0; // C01
        let mut temp: digit_t = 0; // C10
        let mut albl: digit_t = 0;
        let mut albh: digit_t = 0;
        let mut ahbl: digit_t = 0;
        let mut ahbh: digit_t = 0;
        let mut res1: digit_t = 0;
        let mut res2: digit_t = 0;
        let mut res3: digit_t = 0;
        let mut carry: digit_t = 0;
        let mut mask_low: digit_t = -(1 as libc::c_int) as digit_t
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        let mut mask_high: digit_t = (-(1 as libc::c_int) as digit_t)
            << (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        al = a & mask_low;
        ah = a
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        bl = b & mask_low;
        bh = b
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        albl = al.wrapping_mul(bl);
        albh = al.wrapping_mul(bh);
        ahbl = ah.wrapping_mul(bl);
        ahbh = ah.wrapping_mul(bh);
        *c.offset(0 as libc::c_int as isize) = albl & mask_low;
        res1 = albl
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        res2 = ahbl & mask_low;
        res3 = albh & mask_low;
        temp = res1.wrapping_add(res2).wrapping_add(res3);
        carry = temp
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        let ref mut fresh0 = *c.offset(0 as libc::c_int as isize);
        *fresh0 ^= temp
            << (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        res1 = ahbl
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        res2 = albh
            >> (::std::mem::size_of::<digit_t>() as libc::c_ulong)
                .wrapping_mul(4 as libc::c_int as libc::c_ulong);
        res3 = ahbh & mask_low;
        temp = res1
            .wrapping_add(res2)
            .wrapping_add(res3)
            .wrapping_add(carry);
        *c.offset(1 as libc::c_int as isize) = temp & mask_low;
        carry = temp & mask_high;
        let ref mut fresh1 = *c.offset(1 as libc::c_int as isize);
        *fresh1 ^= (ahbh & mask_high).wrapping_add(carry);
        // C11
    }

    pub unsafe extern "C" fn mp_mul(
        mut a: *const digit_t,
        mut b: *const digit_t,
        mut c: *mut digit_t,
        nwords: libc::c_uint,
    ) {
        // Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.
        let mut i: libc::c_uint = 0;
        let mut j: libc::c_uint = 0;
        let mut t: digit_t = 0 as libc::c_int as digit_t;
        let mut u: digit_t = 0 as libc::c_int as digit_t;
        let mut v: digit_t = 0 as libc::c_int as digit_t;
        let mut UV: [digit_t; 2] = [0; 2];
        let mut carry: libc::c_uint = 0 as libc::c_int as libc::c_uint;
        i = 0 as libc::c_int as libc::c_uint;
        while i < nwords {
            j = 0 as libc::c_int as libc::c_uint;
            while j <= i {
                let mut tempReg: uint128_t = (*a.offset(j as isize) as uint128_t)
                    .wrapping_mul(*b.offset(i.wrapping_sub(j) as isize) as uint128_t);
                *UV.as_mut_ptr().offset(1 as libc::c_int as isize) =
                    (tempReg >> 64 as libc::c_int) as digit_t;
                UV[0 as libc::c_int as usize] = tempReg as digit_t;
                let mut tempReg_0: uint128_t = (UV[0 as libc::c_int as usize] as uint128_t)
                    .wrapping_add(v as uint128_t)
                    .wrapping_add(0 as libc::c_int as uint128_t);
                carry = (tempReg_0 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                v = tempReg_0 as digit_t;
                let mut tempReg_1: uint128_t = (UV[1 as libc::c_int as usize] as uint128_t)
                    .wrapping_add(u as uint128_t)
                    .wrapping_add(carry as uint128_t);
                carry = (tempReg_1 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                u = tempReg_1 as digit_t;
                t = (t as libc::c_ulong).wrapping_add(carry as libc::c_ulong) as digit_t as digit_t;
                j = j.wrapping_add(1)
            }
            *c.offset(i as isize) = v;
            v = u;
            u = t;
            t = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
        i = nwords;
        while i
            < (2 as libc::c_int as libc::c_uint)
                .wrapping_mul(nwords)
                .wrapping_sub(1 as libc::c_int as libc::c_uint)
        {
            j = i
                .wrapping_sub(nwords)
                .wrapping_add(1 as libc::c_int as libc::c_uint);
            while j < nwords {
                let mut tempReg_2: uint128_t = (*a.offset(j as isize) as uint128_t)
                    .wrapping_mul(*b.offset(i.wrapping_sub(j) as isize) as uint128_t);
                *UV.as_mut_ptr().offset(1 as libc::c_int as isize) =
                    (tempReg_2 >> 64 as libc::c_int) as digit_t;
                UV[0 as libc::c_int as usize] = tempReg_2 as digit_t;
                let mut tempReg_3: uint128_t = (UV[0 as libc::c_int as usize] as uint128_t)
                    .wrapping_add(v as uint128_t)
                    .wrapping_add(0 as libc::c_int as uint128_t);
                carry = (tempReg_3 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                v = tempReg_3 as digit_t;
                let mut tempReg_4: uint128_t = (UV[1 as libc::c_int as usize] as uint128_t)
                    .wrapping_add(u as uint128_t)
                    .wrapping_add(carry as uint128_t);
                carry = (tempReg_4 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                u = tempReg_4 as digit_t;
                t = (t as libc::c_ulong).wrapping_add(carry as libc::c_ulong) as digit_t as digit_t;
                j = j.wrapping_add(1)
            }
            *c.offset(i as isize) = v;
            v = u;
            u = t;
            t = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
        *c.offset(
            (2 as libc::c_int as libc::c_uint)
                .wrapping_mul(nwords)
                .wrapping_sub(1 as libc::c_int as libc::c_uint) as isize,
        ) = v;
    }

    pub unsafe extern "C" fn rdc_mont(mut ma: *const digit_t, mut mc: *mut digit_t) {
        // Efficient Montgomery reduction using comba and exploiting the special form of the prime p751.
        // mc = ma*R^-1 mod p751x2, where R = 2^768.
        // If ma < 2^768*p751, the output mc is in the range [0, 2*p751-1].
        // ma is assumed to be in Montgomery representation.
        let mut i: libc::c_uint = 0;
        let mut j: libc::c_uint = 0;
        let mut carry: libc::c_uint = 0;
        let mut count: libc::c_uint = 5 as libc::c_int as libc::c_uint;
        let mut UV: [digit_t; 2] = [0; 2];
        let mut t: digit_t = 0 as libc::c_int as digit_t;
        let mut u: digit_t = 0 as libc::c_int as digit_t;
        let mut v: digit_t = 0 as libc::c_int as digit_t;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            *mc.offset(i as isize) = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            j = 0 as libc::c_int as libc::c_uint;
            while j < i {
                if j < i
                    .wrapping_sub(5 as libc::c_int as libc::c_uint)
                    .wrapping_add(1 as libc::c_int as libc::c_uint)
                {
                    let mut tempReg: uint128_t = (*mc.offset(j as isize) as uint128_t)
                        .wrapping_mul(
                            *(crate::p751p1.as_ptr() as *mut digit_t)
                                .offset(i.wrapping_sub(j) as isize)
                                as uint128_t,
                        );
                    *UV.as_mut_ptr().offset(1 as libc::c_int as isize) =
                        (tempReg >> 64 as libc::c_int) as digit_t;
                    UV[0 as libc::c_int as usize] = tempReg as digit_t;
                    let mut tempReg_0: uint128_t = (UV[0 as libc::c_int as usize] as uint128_t)
                        .wrapping_add(v as uint128_t)
                        .wrapping_add(0 as libc::c_int as uint128_t);
                    carry = (tempReg_0 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                    v = tempReg_0 as digit_t;
                    let mut tempReg_1: uint128_t = (UV[1 as libc::c_int as usize] as uint128_t)
                        .wrapping_add(u as uint128_t)
                        .wrapping_add(carry as uint128_t);
                    carry = (tempReg_1 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                    u = tempReg_1 as digit_t;
                    t = (t as libc::c_ulong).wrapping_add(carry as libc::c_ulong) as digit_t
                        as digit_t
                }
                j = j.wrapping_add(1)
            }
            let mut tempReg_2: uint128_t = (v as uint128_t)
                .wrapping_add(*ma.offset(i as isize) as uint128_t)
                .wrapping_add(0 as libc::c_int as uint128_t);
            carry = (tempReg_2 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            v = tempReg_2 as digit_t;
            let mut tempReg_3: uint128_t = (u as uint128_t)
                .wrapping_add(0 as libc::c_int as uint128_t)
                .wrapping_add(carry as uint128_t);
            carry = (tempReg_3 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            u = tempReg_3 as digit_t;
            t = (t as libc::c_ulong).wrapping_add(carry as libc::c_ulong) as digit_t as digit_t;
            *mc.offset(i as isize) = v;
            v = u;
            u = t;
            t = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
        i = 12 as libc::c_int as libc::c_uint;
        while i < (2 as libc::c_int * 12 as libc::c_int - 1 as libc::c_int) as libc::c_uint {
            if count > 0 as libc::c_int as libc::c_uint {
                count = count.wrapping_sub(1 as libc::c_int as libc::c_uint)
            }
            j = i
                .wrapping_sub(12 as libc::c_int as libc::c_uint)
                .wrapping_add(1 as libc::c_int as libc::c_uint);
            while j < 12 as libc::c_int as libc::c_uint {
                if j < (12 as libc::c_int as libc::c_uint).wrapping_sub(count) {
                    let mut tempReg_4: uint128_t = (*mc.offset(j as isize) as uint128_t)
                        .wrapping_mul(
                            *(crate::p751p1.as_ptr() as *mut digit_t)
                                .offset(i.wrapping_sub(j) as isize)
                                as uint128_t,
                        );
                    *UV.as_mut_ptr().offset(1 as libc::c_int as isize) =
                        (tempReg_4 >> 64 as libc::c_int) as digit_t;
                    UV[0 as libc::c_int as usize] = tempReg_4 as digit_t;
                    let mut tempReg_5: uint128_t = (UV[0 as libc::c_int as usize] as uint128_t)
                        .wrapping_add(v as uint128_t)
                        .wrapping_add(0 as libc::c_int as uint128_t);
                    carry = (tempReg_5 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                    v = tempReg_5 as digit_t;
                    let mut tempReg_6: uint128_t = (UV[1 as libc::c_int as usize] as uint128_t)
                        .wrapping_add(u as uint128_t)
                        .wrapping_add(carry as uint128_t);
                    carry = (tempReg_6 >> 64 as libc::c_int) as digit_t as libc::c_uint;
                    u = tempReg_6 as digit_t;
                    t = (t as libc::c_ulong).wrapping_add(carry as libc::c_ulong) as digit_t
                        as digit_t
                }
                j = j.wrapping_add(1)
            }
            let mut tempReg_7: uint128_t = (v as uint128_t)
                .wrapping_add(*ma.offset(i as isize) as uint128_t)
                .wrapping_add(0 as libc::c_int as uint128_t);
            carry = (tempReg_7 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            v = tempReg_7 as digit_t;
            let mut tempReg_8: uint128_t = (u as uint128_t)
                .wrapping_add(0 as libc::c_int as uint128_t)
                .wrapping_add(carry as uint128_t);
            carry = (tempReg_8 >> 64 as libc::c_int) as digit_t as libc::c_uint;
            u = tempReg_8 as digit_t;
            t = (t as libc::c_ulong).wrapping_add(carry as libc::c_ulong) as digit_t as digit_t;
            *mc.offset(i.wrapping_sub(12 as libc::c_int as libc::c_uint) as isize) = v;
            v = u;
            u = t;
            t = 0 as libc::c_int as digit_t;
            i = i.wrapping_add(1)
        }
        let mut tempReg_9: uint128_t = (v as uint128_t)
            .wrapping_add(
                *ma.offset((2 as libc::c_int * 12 as libc::c_int - 1 as libc::c_int) as isize)
                    as uint128_t,
            )
            .wrapping_add(0 as libc::c_int as uint128_t);
        carry = (tempReg_9 >> 64 as libc::c_int) as digit_t as libc::c_uint;
        v = tempReg_9 as digit_t;
        *mc.offset((12 as libc::c_int - 1 as libc::c_int) as isize) = v;
    }
}
