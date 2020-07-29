use crate::fpx::{
    copy_words, fp2add, fp2copy, fp2correction, fp2div2, fp2inv_mont, fp2mul_mont, fp2sqr_mont,
    fp2sub,
};
use crate::{f2elm_t, PointProj};

pub fn xDBL(p: &PointProj, q: &mut PointProj, a24plus: &f2elm_t, c24: &f2elm_t) {
    let (mut t0, mut t1) = ([[0u64; 12]; 2], [[0u64; 12]; 2]);

    fp2sub(&p.x, &p.z, &mut t0);
    fp2add(&p.x, &p.z, &mut t1);
    let mut temp = t0;
    fp2sqr_mont(&temp, &mut t0);
    temp = t1;
    fp2sqr_mont(&temp, &mut t1);
    fp2mul_mont(c24, &t0, &mut q.z);
    fp2mul_mont(&t1, &q.z, &mut q.x);
    temp = t1;
    fp2sub(&temp, &t0, &mut t1);
    fp2mul_mont(a24plus, &t1, &mut t0);
    temp = q.z;
    fp2add(&temp, &t0, &mut q.z);
    temp = q.z;
    fp2mul_mont(&temp, &t1, &mut q.z);
}

pub fn xDBLe(p: &PointProj, q: &mut PointProj, a24plus: &f2elm_t, c24: &f2elm_t, e: usize) {
    copy_words(&p.x[0], &mut q.x[0], 2 * 2 * crate::NWORDS_FIELD);

    for _ in 0..e {
        let temp = *q;
        xDBL(&temp, q, a24plus, c24);
    }
}

pub fn get_4_isog(p: &PointProj, a24plus: &mut f2elm_t, c24: &mut f2elm_t, coeff: &mut [f2elm_t]) {
    fp2sub(&p.x, &p.z, &mut coeff[1]);
    fp2add(&p.x, &p.z, &mut coeff[2]);
    fp2sqr_mont(&p.z, &mut coeff[0]);
    let mut temp = coeff[0];
    fp2add(&temp, &temp, &mut coeff[0]);
    fp2sqr_mont(&coeff[0], c24);
    temp = coeff[0];
    fp2add(&temp, &temp, &mut coeff[0]);
    fp2sqr_mont(&p.x, a24plus);
    temp = *a24plus;
    fp2add(&temp, &temp, a24plus);
    temp = *a24plus;
    fp2sqr_mont(&temp, a24plus);
}

pub fn eval_4_isog(p: &mut PointProj, coeff: &mut [f2elm_t]) {
    let (mut t0, mut t1) = ([[0u64; 12]; 2], [[0u64; 12]; 2]);

    fp2add(&p.x, &p.z, &mut t0);
    fp2sub(&p.x, &p.z, &mut t1);
    fp2mul_mont(&t0, &coeff[1], &mut p.x);
    fp2mul_mont(&t1, &coeff[2], &mut p.z);
    let mut temp0 = t0;
    fp2mul_mont(&temp0, &t1, &mut t0);
    temp0 = t0;
    fp2mul_mont(&temp0, &coeff[0], &mut t0);
    fp2add(&p.x, &p.z, &mut t1);
    let temp1 = *p;
    fp2sub(&temp1.x, &temp1.z, &mut p.z);
    temp0 = t1;
    fp2sqr_mont(&temp0, &mut t1);
    temp0 = p.z;
    fp2sqr_mont(&temp0, &mut p.z);
    fp2add(&t1, &t0, &mut p.x);
    temp0 = t0;
    fp2sub(&p.z, &temp0, &mut t0);
    temp0 = p.x;
    fp2mul_mont(&temp0, &t1, &mut p.x);
    temp0 = p.z;
    fp2mul_mont(&temp0, &t0, &mut p.z);
}

pub fn xTPL(p: &PointProj, q: &mut PointProj, a24minus: &f2elm_t, a24plus: &f2elm_t) {
    let (mut t0, mut t1, mut t2, mut t3, mut t4, mut t5, mut t6) = (
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
    );

    fp2sub(&p.x, &p.z, &mut t0);
    fp2sqr_mont(&t0, &mut t2);
    fp2add(&p.x, &p.z, &mut t1);
    fp2sqr_mont(&t1, &mut t3);
    fp2add(&t0, &t1, &mut t4);
    let mut temp = t0;
    fp2sub(&t1, &temp, &mut t0);
    fp2sqr_mont(&t4, &mut t1);
    temp = t1;
    fp2sub(&temp, &t3, &mut t1);
    temp = t1;
    fp2sub(&temp, &t2, &mut t1);
    fp2mul_mont(&t3, a24plus, &mut t5);
    temp = t3;
    fp2mul_mont(&temp, &t5, &mut t3);
    fp2mul_mont(a24minus, &t2, &mut t6);
    temp = t2;
    fp2mul_mont(&temp, &t6, &mut t2);
    temp = t3;
    fp2sub(&t2, &temp, &mut t3);
    fp2sub(&t5, &t6, &mut t2);
    temp = t1;
    fp2mul_mont(&temp, &t2, &mut t1);
    fp2add(&t3, &t1, &mut t2);
    temp = t2;
    fp2sqr_mont(&temp, &mut t2);
    fp2mul_mont(&t4, &t2, &mut q.x);
    temp = t1;
    fp2sub(&t3, &temp, &mut t1);
    temp = t1;
    fp2sqr_mont(&temp, &mut t1);
    fp2mul_mont(&t0, &t1, &mut q.z);
}

pub fn xTPLe(p: &PointProj, q: &mut PointProj, a24minus: &f2elm_t, a24plus: &f2elm_t, e: usize) {
    copy_words(&p.x[0], &mut q.x[0], 2 * 2 * crate::NWORDS_FIELD);

    for _ in 0..e {
        let temp = *q;
        xTPL(&temp, q, a24minus, a24plus);
    }
}

pub fn get_3_isog(
    p: &PointProj,
    a24minus: &mut f2elm_t,
    a24plus: &mut f2elm_t,
    coeff: &mut [f2elm_t],
) {
    let (mut t0, mut t1, mut t2, mut t3, mut t4) = (
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
    );

    fp2sub(&p.x, &p.z, &mut coeff[0]);
    fp2sqr_mont(&coeff[0], &mut t0);
    fp2add(&p.x, &p.z, &mut coeff[1]);
    fp2sqr_mont(&coeff[1], &mut t1);
    fp2add(&t0, &t1, &mut t2);
    fp2add(&coeff[0], &coeff[1], &mut t3);
    let mut temp = t3;
    fp2sqr_mont(&temp, &mut t3);
    temp = t3;
    fp2sub(&temp, &t2, &mut t3);
    fp2add(&t1, &t3, &mut t2);
    temp = t3;
    fp2add(&temp, &t0, &mut t3);
    fp2add(&t0, &t3, &mut t4);
    temp = t4;
    fp2add(&temp, &temp, &mut t4);
    temp = t4;
    fp2add(&t1, &temp, &mut t4);
    fp2mul_mont(&t2, &t4, a24minus);
    fp2add(&t1, &t2, &mut t4);
    temp = t4;
    fp2add(&temp, &temp, &mut t4);
    temp = t4;
    fp2add(&t0, &temp, &mut t4);
    fp2mul_mont(&t3, &t4, a24plus);
}

pub fn eval_3_isog(q: &mut PointProj, coeff: &[f2elm_t]) {
    let (mut t0, mut t1, mut t2) = ([[0u64; 12]; 2], [[0u64; 12]; 2], [[0u64; 12]; 2]);

    fp2add(&q.x, &q.z, &mut t0);
    fp2sub(&q.x, &q.z, &mut t1);
    let mut temp = t0;
    fp2mul_mont(&temp, &coeff[0], &mut t0);
    temp = t1;
    fp2mul_mont(&temp, &coeff[1], &mut t1);
    fp2add(&t0, &t1, &mut t2);
    temp = t0;
    fp2sub(&t1, &temp, &mut t0);
    temp = t2;
    fp2sqr_mont(&temp, &mut t2);
    temp = t0;
    fp2sqr_mont(&temp, &mut t0);
    temp = q.x;
    fp2mul_mont(&temp, &t2, &mut q.x);
    temp = q.z;
    fp2mul_mont(&temp, &t0, &mut q.z);
}

pub fn inv_3_way(z1: &mut f2elm_t, z2: &mut f2elm_t, z3: &mut f2elm_t) {
    let (mut t0, mut t1, mut t2, mut t3) = (
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
        [[0u64; 12]; 2],
    );

    fp2mul_mont(z1, z2, &mut t0);
    fp2mul_mont(z3, &t0, &mut t1);
    fp2inv_mont(&mut t1);
    fp2mul_mont(z3, &t1, &mut t2);
    fp2mul_mont(&t2, z2, &mut t3);
    fp2mul_mont(&t2, z1, z2);
    fp2mul_mont(&t0, &t1, z3);
    fp2copy(&t3, z1);
}

pub fn get_a(xP: &f2elm_t, xQ: &f2elm_t, xR: &f2elm_t, a: &mut f2elm_t) {
    let (mut t0, mut t1, mut one) = ([[0u64; 12]; 2], [[0u64; 12]; 2], [[0u64; 12]; 2]);

    one[0] = crate::Montgomery_one;
    fp2add(xP, xQ, &mut t1);
    fp2mul_mont(xP, xQ, &mut t0);
    fp2mul_mont(xR, &t1, a);
    let mut temp = *a;
    fp2add(&t0, &temp, a);
    temp = t0;
    fp2mul_mont(&temp, xR, &mut t0);
    temp = *a;
    fp2sub(&temp, &one, a);
    temp = t0;
    fp2add(&temp, &temp, &mut t0);
    temp = t1;
    fp2add(&temp, xR, &mut t1);
    temp = t0;
    fp2add(&temp, &temp, &mut t0);
    temp = *a;
    fp2sqr_mont(&temp, a);
    fp2inv_mont(&mut t0);
    temp = *a;
    fp2mul_mont(&temp, &t0, a);
    temp = *a;
    fp2sub(&temp, &t1, a);
}

pub fn j_inv(a: &f2elm_t, c: &f2elm_t, jinv: &mut f2elm_t) {
    let (mut t0, mut t1) = ([[0u64; 12]; 2], [[0u64; 12]; 2]);

    fp2sqr_mont(a, jinv);
    fp2sqr_mont(c, &mut t1);
    fp2add(&t1, &t1, &mut t0);
    let mut temp = t0;
    fp2sub(jinv, &temp, &mut t0);
    temp = t0;
    fp2sub(&temp, &t1, &mut t0);
    fp2sub(&t0, &t1, jinv);
    temp = t1;
    fp2sqr_mont(&temp, &mut t1);
    temp = *jinv;
    fp2mul_mont(&temp, &t1, jinv);
    temp = t0;
    fp2add(&temp, &temp, &mut t0);
    temp = t0;
    fp2add(&temp, &temp, &mut t0);
    fp2sqr_mont(&t0, &mut t1);
    temp = t0;
    fp2mul_mont(&temp, &t1, &mut t0);
    temp = t0;
    fp2add(&temp, &temp, &mut t0);
    temp = t0;
    fp2add(&temp, &temp, &mut t0);
    fp2inv_mont(jinv);
    temp = *jinv;
    fp2mul_mont(&temp, &t0, jinv);
}

pub fn xDBLaDD(p: &mut PointProj, q: &mut PointProj, xPQ: &f2elm_t, a24: &f2elm_t) {
    let (mut t0, mut t1, mut t2) = ([[0u64; 12]; 2], [[0u64; 12]; 2], [[0u64; 12]; 2]);

    fp2add(&p.x, &p.z, &mut t0);
    fp2sub(&p.x, &p.z, &mut t1);
    fp2sqr_mont(&t0, &mut p.x);
    fp2sub(&q.x, &q.z, &mut t2);
    fp2correction(&mut t2);
    let mut temp = q.x;
    fp2add(&temp, &q.z, &mut q.x);
    temp = t0;
    fp2mul_mont(&temp, &t2, &mut t0);
    fp2sqr_mont(&t1, &mut p.z);
    temp = t1;
    fp2mul_mont(&temp, &q.x, &mut t1);
    fp2sub(&p.x, &p.z, &mut t2);
    temp = p.x;
    fp2mul_mont(&temp, &p.z, &mut p.x);
    fp2mul_mont(&t2, a24, &mut q.x);
    fp2sub(&t0, &t1, &mut q.z);
    temp = p.z;
    fp2add(&q.x, &temp, &mut p.z);
    fp2add(&t0, &t1, &mut q.x);
    temp = p.z;
    fp2mul_mont(&temp, &t2, &mut p.z);
    temp = q.z;
    fp2sqr_mont(&temp, &mut q.z);
    temp = q.x;
    fp2sqr_mont(&temp, &mut q.x);
    temp = q.z;
    fp2mul_mont(&temp, xPQ, &mut q.z);
}

pub fn swap_points(p: &mut PointProj, q: &mut PointProj, option: u64) {
    for i in 0..12 as usize {
        let mut temp = option & (p.x[0][i] ^ q.x[0][i]);
        p.x[0][i] = temp ^ p.x[0][i];
        q.x[0][i] = temp ^ q.x[0][i];
        temp = option & (p.z[0][i] ^ q.z[0][i]);
        p.z[0][i] = temp ^ p.z[0][i];
        q.z[0][i] = temp ^ q.z[0][i];
        temp = option & (p.x[1][i] ^ q.x[1][i]);
        p.x[1][i] = temp ^ p.x[1][i];
        q.x[1][i] = temp ^ q.x[1][i];
        temp = option & (p.z[1][i] ^ q.z[1][i]);
        p.z[1][i] = temp ^ p.z[1][i];
        q.z[1][i] = temp ^ q.z[1][i];
    }
}

pub fn ladder3pt(
    xP: &f2elm_t,
    xQ: &f2elm_t,
    xPQ: &f2elm_t,
    m: &[u64],
    alice_or_bob: u32,
    r: &mut PointProj,
    a: &f2elm_t,
) {
    let mut R0 = PointProj::new();
    let mut R2 = PointProj::new();
    let mut a24 = [[0u64; 12]; 2];
    let mut nbits: i32 = 0;
    let mut prevbit: i32 = 0;

    if alice_or_bob == 0 {
        nbits = 372;
    }

    a24[0] = crate::Montgomery_one;

    let mut temp = a24;
    fp2add(&temp, &temp, &mut a24);
    temp = a24;
    fp2add(a, &temp, &mut a24);
    temp = a24;
    fp2div2(&temp, &mut a24);
    temp = a24;
    fp2div2(&temp, &mut a24);

    R0.x = *xQ;
    R0.z[0] = crate::Montgomery_one;
    R2.x = *xPQ;
    R2.z[0] = crate::Montgomery_one;
    r.x = *xP;
    r.z[0] = crate::Montgomery_one;
    r.z[1].iter_mut().for_each(|x| *x = 0);

    for i in 0..nbits {
        let bit = (m[(i >> crate::LOG2RADIX as libc::c_int) as usize]
            >> (i & crate::RADIX as i32 - 1)
            & 1 as u64) as i32;

        prevbit = bit;

        swap_points(
            r,
            &mut R2,
            (0 as u64).wrapping_sub((0 as u64).wrapping_sub((bit ^ prevbit) as u64)),
        );
        xDBLaDD(&mut R0, &mut R2, &r.x, &a24);
        let temp = R2.x;
        fp2mul_mont(&temp, &r.z, &mut R2.x);
    }

    swap_points(r, &mut R2, (0 as u64).wrapping_sub((0 ^ prevbit) as u64));
}
