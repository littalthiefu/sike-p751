use crate::fpx::{
    copy_words, copy_words_x, fp2add, fp2addx, fp2copy, fp2correction, fp2div2, fp2inv_mont,
    fp2mul_mont, fp2mul_montx, fp2sqr_mont, fp2sqr_montx, fp2sub, fp2subx,
};
use crate::{f2elm_t, PointProj, PointProjRaw};

pub fn xDBL(p: *const PointProj, q: &mut PointProj, a24plus: &f2elm_t, c24: &f2elm_t) {
    let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = X1-Z1
    let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = X1+Z1

    unsafe {
        fp2sub(&(*p).x, &(*p).z, &mut t0); // t0 = (X1-Z1)^2
        fp2add(&(*p).x, &(*p).z, &mut t1); // t1 = (X1+Z1)^2
        fp2sqr_montx(t0.as_ptr(), &mut t0); // Z2 = C24*(X1-Z1)^2
        fp2sqr_montx(t1.as_ptr(), &mut t1); // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
        fp2mul_mont(c24, &t0, &mut q.z); // t1 = (X1+Z1)^2-(X1-Z1)^2
        fp2mul_montx(t1.as_ptr(), q.z.as_ptr(), &mut q.x); // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
        fp2subx(t1.as_ptr(), t0.as_ptr(), &mut t1); // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
        fp2mul_mont(a24plus, &t1, &mut t0);
        fp2addx(q.z.as_ptr(), t0.as_ptr(), &mut q.z);
        fp2mul_montx(q.z.as_ptr(), t1.as_ptr(), &mut q.z);
    }
}

pub fn xDBLe(p: *const PointProj, q: &mut PointProj, a24plus: &f2elm_t, c24: &f2elm_t, e: usize) {
    let p0 = PointProjRaw::from(unsafe { *p });
    let mut q0 = PointProjRaw::from(unsafe { *q });

    copy_words(
        &p0.0,
        &mut q0.0,
        (2 as libc::c_int * 2 as libc::c_int * 12 as libc::c_int) as usize,
    );

    let mut q0 = [PointProj::from(q0.0)];

    for _ in 0..e {
        xDBL(q0.as_ptr(), &mut q0[0], a24plus, c24);
    }

    *q = q0[0];
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
    m: *const u64,
    alice_or_bob: usize,
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
        let bit = unsafe {
            (*m.offset((i >> crate::LOG2RADIX as libc::c_int) as isize)
                >> (i & crate::RADIX as i32 - 1)
                & 1 as u64) as i32
        };

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

#[cfg(test)]
pub mod tests {
    use crate::fpx::tests::*;
    use crate::*;

    #[no_mangle]
    pub unsafe extern "C" fn xDBL(
        mut P: *const point_proj,
        mut Q: *mut point_proj,
        mut A24plus: *const felm_t,
        mut C24: *const felm_t,
    ) {
        // Doubling of a Montgomery point in projective coordinates (X:Z).
        // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.
        // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
        let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = X1-Z1
        let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = X1+Z1
        fp2sub((*P).X.as_ptr(), (*P).Z.as_ptr(), t0.as_mut_ptr()); // t0 = (X1-Z1)^2
        fp2add((*P).X.as_ptr(), (*P).Z.as_ptr(), t1.as_mut_ptr()); // t1 = (X1+Z1)^2
        fp2sqr_mont(t0.as_mut_ptr() as *const felm_t, t0.as_mut_ptr()); // Z2 = C24*(X1-Z1)^2
        fp2sqr_mont(t1.as_mut_ptr() as *const felm_t, t1.as_mut_ptr()); // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
        fp2mul_mont(C24, t0.as_mut_ptr() as *const felm_t, (*Q).Z.as_mut_ptr()); // t1 = (X1+Z1)^2-(X1-Z1)^2
        fp2mul_mont(
            t1.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr() as *const felm_t,
            (*Q).X.as_mut_ptr(),
        ); // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
        fp2sub(
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
        fp2mul_mont(A24plus, t1.as_mut_ptr() as *const felm_t, t0.as_mut_ptr());
        fp2add(
            (*Q).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr(),
        );
        fp2mul_mont(
            (*Q).Z.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr(),
        );
        // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
    }
    #[no_mangle]
    pub unsafe extern "C" fn xDBLe(
        mut P: *const point_proj,
        mut Q: *mut point_proj,
        mut A24plus: *const felm_t,
        mut C24: *const felm_t,
        e: libc::c_int,
    ) {
        // Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
        // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A+2C and 4C.
        // Output: projective Montgomery x-coordinates Q <- (2^e)*P.
        let mut i: libc::c_int = 0;
        copy_words(
            P as *mut digit_t,
            Q as *mut digit_t,
            (2 as libc::c_int * 2 as libc::c_int * 12 as libc::c_int) as libc::c_uint,
        );
        i = 0 as libc::c_int;
        while i < e {
            xDBL(Q as *const point_proj, Q, A24plus, C24);
            i += 1
        }
    }
    #[no_mangle]
    pub unsafe extern "C" fn get_4_isog(
        mut P: *const point_proj,
        mut A24plus: *mut felm_t,
        mut C24: *mut felm_t,
        mut coeff: *mut f2elm_t,
    ) {
        // Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
        // Input:  projective point of order four P = (X4:Z4).
        // Output: the 4-isogenous Montgomery curve with projective coefficients A+2C/4C and the 3 coefficients
        //         that are used to evaluate the isogeny at a point in eval_4_isog().
        fp2sub(
            (*P).X.as_ptr(),
            (*P).Z.as_ptr(),
            (*coeff.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        ); // coeff[1] = X4-Z4
        fp2add(
            (*P).X.as_ptr(),
            (*P).Z.as_ptr(),
            (*coeff.offset(2 as libc::c_int as isize)).as_mut_ptr(),
        ); // coeff[2] = X4+Z4
        fp2sqr_mont(
            (*P).Z.as_ptr(),
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        ); // coeff[0] = Z4^2
        fp2add(
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        ); // coeff[0] = 2*Z4^2
        fp2sqr_mont(
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            C24,
        ); // C24 = 4*Z4^4
        fp2add(
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        ); // coeff[0] = 4*Z4^2
        fp2sqr_mont((*P).X.as_ptr(), A24plus); // A24plus = X4^2
        fp2add(A24plus as *const felm_t, A24plus as *const felm_t, A24plus); // A24plus = 2*X4^2
        fp2sqr_mont(A24plus as *const felm_t, A24plus);
        // A24plus = 4*X4^4
    }
    #[no_mangle]
    pub unsafe extern "C" fn eval_4_isog(mut P: *mut point_proj, mut coeff: *mut f2elm_t) {
        // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined
        // by the 3 coefficients in coeff (computed in the function get_4_isog()).
        // Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
        // Output: the projective point P = phi(P) = (X:Z) in the codomain.
        let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = X+Z
        let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = X-Z
        fp2add(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // X = (X+Z)*coeff[1]
        fp2sub(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // Z = (X-Z)*coeff[2]
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            (*coeff.offset(1 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*P).X.as_mut_ptr(),
        ); // t0 = (X+Z)*(X-Z)
        fp2mul_mont(
            t1.as_mut_ptr() as *const felm_t,
            (*coeff.offset(2 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr(),
        ); // t0 = coeff[0]*(X+Z)*(X-Z)
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
        fp2add(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
        fp2sub(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr(),
        ); // Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
        fp2sqr_mont(t1.as_mut_ptr() as *const felm_t, t1.as_mut_ptr()); // X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
        fp2sqr_mont((*P).Z.as_mut_ptr() as *const felm_t, (*P).Z.as_mut_ptr()); // t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z)
        fp2add(
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            (*P).X.as_mut_ptr(),
        ); // Xfinal
        fp2sub(
            (*P).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        );
        fp2mul_mont(
            (*P).X.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            (*P).X.as_mut_ptr(),
        );
        fp2mul_mont(
            (*P).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr(),
        );
        // Zfinal
    }
    #[no_mangle]
    pub unsafe extern "C" fn xTPL(
        mut P: *const point_proj,
        mut Q: *mut point_proj,
        mut A24minus: *const felm_t,
        mut A24plus: *const felm_t,
    ) {
        // Tripling of a Montgomery point in projective coordinates (X:Z).
        // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
        // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
        let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = X-Z
        let mut t1: f2elm_t = [[0; 12]; 2]; // t2 = (X-Z)^2
        let mut t2: f2elm_t = [[0; 12]; 2]; // t1 = X+Z
        let mut t3: f2elm_t = [[0; 12]; 2]; // t3 = (X+Z)^2
        let mut t4: f2elm_t = [[0; 12]; 2]; // t4 = 2*X
        let mut t5: f2elm_t = [[0; 12]; 2]; // t0 = 2*Z
        let mut t6: f2elm_t = [[0; 12]; 2]; // t1 = 4*X^2
        fp2sub((*P).X.as_ptr(), (*P).Z.as_ptr(), t0.as_mut_ptr()); // t1 = 4*X^2 - (X+Z)^2
        fp2sqr_mont(t0.as_mut_ptr() as *const felm_t, t2.as_mut_ptr()); // t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
        fp2add((*P).X.as_ptr(), (*P).Z.as_ptr(), t1.as_mut_ptr()); // t5 = A24plus*(X+Z)^2
        fp2sqr_mont(t1.as_mut_ptr() as *const felm_t, t3.as_mut_ptr()); // t3 = A24plus*(X+Z)^3
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        ); // t6 = A24minus*(X-Z)^2
        fp2sub(
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t2 = A24minus*(X-Z)^3
        fp2sqr_mont(t4.as_mut_ptr() as *const felm_t, t1.as_mut_ptr()); // t3 = A24minus*(X-Z)^3 - coeff*(X+Z)^3
        fp2sub(
            t1.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
        fp2sub(
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
        fp2mul_mont(t3.as_mut_ptr() as *const felm_t, A24plus, t5.as_mut_ptr()); // t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^3 - coeff*(X+Z)^3
        fp2mul_mont(
            t3.as_mut_ptr() as *const felm_t,
            t5.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr(),
        ); // t2 = t2^2
        fp2mul_mont(A24minus, t2.as_mut_ptr() as *const felm_t, t6.as_mut_ptr()); // X3 = 2*X*t2
        fp2mul_mont(
            t2.as_mut_ptr() as *const felm_t,
            t6.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        ); // t1 = A24minus*(X-Z)^3 - A24plus*(X+Z)^3 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
        fp2sub(
            t2.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr(),
        ); // t1 = t1^2
        fp2sub(
            t5.as_mut_ptr() as *const felm_t,
            t6.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        );
        fp2mul_mont(
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        );
        fp2add(
            t3.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        );
        fp2sqr_mont(t2.as_mut_ptr() as *const felm_t, t2.as_mut_ptr());
        fp2mul_mont(
            t4.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            (*Q).X.as_mut_ptr(),
        );
        fp2sub(
            t3.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        );
        fp2sqr_mont(t1.as_mut_ptr() as *const felm_t, t1.as_mut_ptr());
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr(),
        );
        // Z3 = 2*Z*t1
    }
    #[no_mangle]
    pub unsafe extern "C" fn xTPLe(
        mut P: *const point_proj,
        mut Q: *mut point_proj,
        mut A24minus: *const felm_t,
        mut A24plus: *const felm_t,
        e: libc::c_int,
    ) {
        // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
        // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
        // Output: projective Montgomery x-coordinates Q <- (3^e)*P.
        let mut i: libc::c_int = 0;
        copy_words(
            P as *mut digit_t,
            Q as *mut digit_t,
            (2 as libc::c_int * 2 as libc::c_int * 12 as libc::c_int) as libc::c_uint,
        );
        i = 0 as libc::c_int;
        while i < e {
            xTPL(Q as *const point_proj, Q, A24minus, A24plus);
            i += 1
        }
    }
    #[no_mangle]
    pub unsafe extern "C" fn get_3_isog(
        mut P: *const point_proj,
        mut A24minus: *mut felm_t,
        mut A24plus: *mut felm_t,
        mut coeff: *mut f2elm_t,
    ) {
        // Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
        // Input:  projective point of order three P = (X3:Z3).
        // Output: the 3-isogenous Montgomery curve with projective coefficient A/C.
        let mut t0: f2elm_t = [[0; 12]; 2]; // coeff0 = X-Z
        let mut t1: f2elm_t = [[0; 12]; 2]; // t0 = (X-Z)^2
        let mut t2: f2elm_t = [[0; 12]; 2]; // coeff1 = X+Z
        let mut t3: f2elm_t = [[0; 12]; 2]; // t1 = (X+Z)^2
        let mut t4: f2elm_t = [[0; 12]; 2]; // t2 = (X+Z)^2 + (X-Z)^2
        fp2sub(
            (*P).X.as_ptr(),
            (*P).Z.as_ptr(),
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr(),
        ); // t3 = 2*X
        fp2sqr_mont(
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t3 = 4*X^2
        fp2add(
            (*P).X.as_ptr(),
            (*P).Z.as_ptr(),
            (*coeff.offset(1 as libc::c_int as isize)).as_mut_ptr(),
        ); // t3 = 4*X^2 - (X+Z)^2 - (X-Z)^2
        fp2sqr_mont(
            (*coeff.offset(1 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // t2 = 4*X^2 - (X-Z)^2
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        ); // t3 = 4*X^2 - (X+Z)^2
        fp2add(
            (*coeff.offset(0 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            (*coeff.offset(1 as libc::c_int as isize)).as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr(),
        ); // t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2
        fp2sqr_mont(t3.as_mut_ptr() as *const felm_t, t3.as_mut_ptr()); // t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2)
        fp2sub(
            t3.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr(),
        ); // t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
        fp2add(
            t1.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        ); // A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
        fp2add(
            t3.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr(),
        ); // t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t3.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        ); // t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2)
        fp2add(
            t4.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        ); // t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
        fp2add(
            t1.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        );
        fp2mul_mont(
            t2.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr() as *const felm_t,
            A24minus,
        );
        fp2add(
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        );
        fp2add(
            t4.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        );
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr(),
        );
        fp2mul_mont(
            t3.as_mut_ptr() as *const felm_t,
            t4.as_mut_ptr() as *const felm_t,
            A24plus,
        );
        // A24plus = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
    }
    #[no_mangle]
    pub unsafe extern "C" fn eval_3_isog(mut Q: *mut point_proj, mut coeff: *const f2elm_t) {
        // Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and
        // a point P with 2 coefficients in coeff (computed in the function get_3_isog()).
        // Inputs: projective points P = (X3:Z3) and Q = (X:Z).
        // Output: the projective point Q <- phi(Q) = (X3:Z3).
        let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = X+Z
        let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = X-Z
        let mut t2: f2elm_t = [[0; 12]; 2]; // t0 = coeff0*(X+Z)
        fp2add(
            (*Q).X.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t1 = coeff1*(X-Z)
        fp2sub(
            (*Q).X.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // t2 = coeff0*(X+Z) + coeff1*(X-Z)
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            (*coeff.offset(0 as libc::c_int as isize)).as_ptr(),
            t0.as_mut_ptr(),
        ); // t0 = coeff1*(X-Z) - coeff0*(X+Z)
        fp2mul_mont(
            t1.as_mut_ptr() as *const felm_t,
            (*coeff.offset(1 as libc::c_int as isize)).as_ptr(),
            t1.as_mut_ptr(),
        ); // t2 = [coeff0*(X+Z) + coeff1*(X-Z)]^2
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        ); // t0 = [coeff1*(X-Z) - coeff0*(X+Z)]^2
        fp2sub(
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // X3final = X*[coeff0*(X+Z) + coeff1*(X-Z)]^2
        fp2sqr_mont(t2.as_mut_ptr() as *const felm_t, t2.as_mut_ptr());
        fp2sqr_mont(t0.as_mut_ptr() as *const felm_t, t0.as_mut_ptr());
        fp2mul_mont(
            (*Q).X.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            (*Q).X.as_mut_ptr(),
        );
        fp2mul_mont(
            (*Q).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr(),
        );
        // Z3final = Z*[coeff1*(X-Z) - coeff0*(X+Z)]^2
    }
    #[no_mangle]
    pub unsafe extern "C" fn inv_3_way(
        mut z1: *mut felm_t,
        mut z2: *mut felm_t,
        mut z3: *mut felm_t,
    ) {
        // 3-way simultaneous inversion
        // Input:  z1,z2,z3
        // Output: 1/z1,1/z2,1/z3 (override inputs).
        let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = z1*z2
        let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = z1*z2*z3
        let mut t2: f2elm_t = [[0; 12]; 2]; // t1 = 1/(z1*z2*z3)
        let mut t3: f2elm_t = [[0; 12]; 2]; // t2 = 1/(z1*z2)
        fp2mul_mont(z1 as *const felm_t, z2 as *const felm_t, t0.as_mut_ptr()); // t3 = 1/z1
        fp2mul_mont(
            z3 as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // z2 = 1/z2
        fp2inv_mont(t1.as_mut_ptr()); // z3 = 1/z3
        fp2mul_mont(
            z3 as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        );
        fp2mul_mont(
            t2.as_mut_ptr() as *const felm_t,
            z2 as *const felm_t,
            t3.as_mut_ptr(),
        );
        fp2mul_mont(t2.as_mut_ptr() as *const felm_t, z1 as *const felm_t, z2);
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            z3,
        );
        fp2copy(t3.as_mut_ptr() as *const felm_t, z1);
        // z1 = 1/z1
    }
    #[no_mangle]
    pub unsafe extern "C" fn get_A(
        mut xP: *const felm_t,
        mut xQ: *const felm_t,
        mut xR: *const felm_t,
        mut A: *mut felm_t,
    ) {
        // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
        // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
        // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
        let mut t0: f2elm_t = [[0; 12]; 2]; // t1 = xP+xQ
        let mut t1: f2elm_t = [[0; 12]; 2]; // t0 = xP*xQ
        let mut one: f2elm_t = [
            [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0; 12],
        ]; // A = xR*t1
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            one[0 as libc::c_int as usize].as_mut_ptr(),
        ); // A = A+t0
        fp2add(xP, xQ, t1.as_mut_ptr()); // t0 = t0*xR
        fp2mul_mont(xP, xQ, t0.as_mut_ptr()); // A = A-1
        fp2mul_mont(xR, t1.as_mut_ptr() as *const felm_t, A); // t0 = t0+t0
        fp2add(t0.as_mut_ptr() as *const felm_t, A as *const felm_t, A); // t1 = t1+xR
        fp2mul_mont(t0.as_mut_ptr() as *const felm_t, xR, t0.as_mut_ptr()); // t0 = t0+t0
        fp2sub(A as *const felm_t, one.as_mut_ptr() as *const felm_t, A); // A = A^2
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t0 = 1/t0
        fp2add(t1.as_mut_ptr() as *const felm_t, xR, t1.as_mut_ptr()); // A = A*t0
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        );
        fp2sqr_mont(A as *const felm_t, A);
        fp2inv_mont(t0.as_mut_ptr());
        fp2mul_mont(A as *const felm_t, t0.as_mut_ptr() as *const felm_t, A);
        fp2sub(A as *const felm_t, t1.as_mut_ptr() as *const felm_t, A);
        // Afinal = A-t1
    }
    #[no_mangle]
    pub unsafe extern "C" fn j_inv(
        mut A: *const felm_t,
        mut C: *const felm_t,
        mut jinv: *mut felm_t,
    ) {
        // Computes the j-invariant of a Montgomery curve with projective constant.
        // Input: A,C in GF(p^2).
        // Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
        let mut t0: f2elm_t = [[0; 12]; 2]; // jinv = A^2
        let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = C^2
        fp2sqr_mont(A, jinv); // t0 = t1+t1
        fp2sqr_mont(C, t1.as_mut_ptr()); // t0 = jinv-t0
        fp2add(
            t1.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t0 = t0-t1
        fp2sub(
            jinv as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // jinv = t0-t1
        fp2sub(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t1 = t1^2
        fp2sub(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            jinv,
        ); // jinv = jinv*t1
        fp2sqr_mont(t1.as_mut_ptr() as *const felm_t, t1.as_mut_ptr()); // t0 = t0+t0
        fp2mul_mont(
            jinv as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            jinv,
        ); // t0 = t0+t0
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t1 = t0^2
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t0 = t0*t1
        fp2sqr_mont(t0.as_mut_ptr() as *const felm_t, t1.as_mut_ptr()); // t0 = t0+t0
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t0 = t0+t0
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // jinv = 1/jinv
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        );
        fp2inv_mont(jinv);
        fp2mul_mont(
            jinv as *const felm_t,
            t0.as_mut_ptr() as *const felm_t,
            jinv,
        );
        // jinv = t0*jinv
    }
    #[no_mangle]
    pub unsafe extern "C" fn xDBLADD(
        mut P: *mut point_proj,
        mut Q: *mut point_proj,
        mut xPQ: *const felm_t,
        mut A24: *const felm_t,
    ) {
        // Simultaneous doubling and differential addition.
        // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
        // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP.
        let mut t0: f2elm_t = [[0; 12]; 2]; // t0 = XP+ZP
        let mut t1: f2elm_t = [[0; 12]; 2]; // t1 = XP-ZP
        let mut t2: f2elm_t = [[0; 12]; 2]; // XP = (XP+ZP)^2
        fp2add(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // t2 = XQ-ZQ
        fp2sub(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // XQ = XQ+ZQ
        fp2sqr_mont(t0.as_mut_ptr() as *const felm_t, (*P).X.as_mut_ptr()); // t0 = (XP+ZP)*(XQ-ZQ)
        fp2sub(
            (*Q).X.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        ); // ZP = (XP-ZP)^2
        fp2correction(t2.as_mut_ptr()); // t1 = (XP-ZP)*(XQ+ZQ)
        fp2add(
            (*Q).X.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr() as *const felm_t,
            (*Q).X.as_mut_ptr(),
        ); // t2 = (XP+ZP)^2-(XP-ZP)^2
        fp2mul_mont(
            t0.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            t0.as_mut_ptr(),
        ); // XP = (XP+ZP)^2*(XP-ZP)^2
        fp2sqr_mont(t1.as_mut_ptr() as *const felm_t, (*P).Z.as_mut_ptr()); // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
        fp2mul_mont(
            t1.as_mut_ptr() as *const felm_t,
            (*Q).X.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr(),
        ); // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
        fp2sub(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr(),
        ); // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
        fp2mul_mont(
            (*P).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            (*P).X.as_mut_ptr(),
        ); // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
        fp2mul_mont(t2.as_mut_ptr() as *const felm_t, A24, (*Q).X.as_mut_ptr()); // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
        fp2sub(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            (*Q).Z.as_mut_ptr(),
        ); // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
        fp2add(
            (*Q).X.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr(),
        ); // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
        fp2add(
            t0.as_mut_ptr() as *const felm_t,
            t1.as_mut_ptr() as *const felm_t,
            (*Q).X.as_mut_ptr(),
        );
        fp2mul_mont(
            (*P).Z.as_mut_ptr() as *const felm_t,
            t2.as_mut_ptr() as *const felm_t,
            (*P).Z.as_mut_ptr(),
        );
        fp2sqr_mont((*Q).Z.as_mut_ptr() as *const felm_t, (*Q).Z.as_mut_ptr());
        fp2sqr_mont((*Q).X.as_mut_ptr() as *const felm_t, (*Q).X.as_mut_ptr());
        fp2mul_mont(
            (*Q).Z.as_mut_ptr() as *const felm_t,
            xPQ,
            (*Q).Z.as_mut_ptr(),
        );
        // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    }
    pub unsafe extern "C" fn swap_points(
        mut P: *mut point_proj,
        mut Q: *mut point_proj,
        option: digit_t,
    ) {
        // Swap points.
        // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
        let mut temp: digit_t = 0;
        let mut i: libc::c_uint = 0;
        i = 0 as libc::c_int as libc::c_uint;
        while i < 12 as libc::c_int as libc::c_uint {
            temp = option
                & ((*P).X[0 as libc::c_int as usize][i as usize]
                    ^ (*Q).X[0 as libc::c_int as usize][i as usize]);
            (*P).X[0 as libc::c_int as usize][i as usize] =
                temp ^ (*P).X[0 as libc::c_int as usize][i as usize];
            (*Q).X[0 as libc::c_int as usize][i as usize] =
                temp ^ (*Q).X[0 as libc::c_int as usize][i as usize];
            temp = option
                & ((*P).Z[0 as libc::c_int as usize][i as usize]
                    ^ (*Q).Z[0 as libc::c_int as usize][i as usize]);
            (*P).Z[0 as libc::c_int as usize][i as usize] =
                temp ^ (*P).Z[0 as libc::c_int as usize][i as usize];
            (*Q).Z[0 as libc::c_int as usize][i as usize] =
                temp ^ (*Q).Z[0 as libc::c_int as usize][i as usize];
            temp = option
                & ((*P).X[1 as libc::c_int as usize][i as usize]
                    ^ (*Q).X[1 as libc::c_int as usize][i as usize]);
            (*P).X[1 as libc::c_int as usize][i as usize] =
                temp ^ (*P).X[1 as libc::c_int as usize][i as usize];
            (*Q).X[1 as libc::c_int as usize][i as usize] =
                temp ^ (*Q).X[1 as libc::c_int as usize][i as usize];
            temp = option
                & ((*P).Z[1 as libc::c_int as usize][i as usize]
                    ^ (*Q).Z[1 as libc::c_int as usize][i as usize]);
            (*P).Z[1 as libc::c_int as usize][i as usize] =
                temp ^ (*P).Z[1 as libc::c_int as usize][i as usize];
            (*Q).Z[1 as libc::c_int as usize][i as usize] =
                temp ^ (*Q).Z[1 as libc::c_int as usize][i as usize];
            i = i.wrapping_add(1)
        }
    }
    pub unsafe extern "C" fn LADDER3PT(
        mut xP: *const felm_t,
        mut xQ: *const felm_t,
        mut xPQ: *const felm_t,
        mut m: *const digit_t,
        AliceOrBob: libc::c_uint,
        mut R: *mut point_proj,
        mut A: *const felm_t,
    ) {
        let mut R0: point_proj_t = [{
            let mut init = point_proj {
                X: [
                    [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0; 12],
                ],
                Z: [[0; 12]; 2],
            };
            init
        }];
        let mut R2: point_proj_t = [{
            let mut init = point_proj {
                X: [
                    [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0; 12],
                ],
                Z: [[0; 12]; 2],
            };
            init
        }];
        let mut A24: f2elm_t = [
            [0 as libc::c_int as digit_t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0; 12],
        ];
        let mut mask: digit_t = 0;
        let mut i: libc::c_int = 0;
        let mut nbits: libc::c_int = 0;
        let mut bit: libc::c_int = 0;
        let mut swap: libc::c_int = 0;
        let mut prevbit: libc::c_int = 0 as libc::c_int;
        if AliceOrBob == 0 as libc::c_int as libc::c_uint {
            nbits = 372 as libc::c_int
        } else {
            nbits = 379 as libc::c_int - 1 as libc::c_int
        }
        // Initializing constant
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            A24[0 as libc::c_int as usize].as_mut_ptr(),
        ); // A24 = (A+2)/4
        fp2add(
            A24.as_mut_ptr() as *const felm_t,
            A24.as_mut_ptr() as *const felm_t,
            A24.as_mut_ptr(),
        );
        fp2add(A, A24.as_mut_ptr() as *const felm_t, A24.as_mut_ptr());
        fp2div2(A24.as_mut_ptr() as *const felm_t, A24.as_mut_ptr());
        fp2div2(A24.as_mut_ptr() as *const felm_t, A24.as_mut_ptr());
        // Initializing points
        fp2copy(xQ, (*R0.as_mut_ptr()).X.as_mut_ptr());
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            (*R0.as_mut_ptr()).Z.as_mut_ptr() as *mut digit_t,
        );
        fp2copy(xPQ, (*R2.as_mut_ptr()).X.as_mut_ptr());
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            (*R2.as_mut_ptr()).Z.as_mut_ptr() as *mut digit_t,
        );
        fp2copy(xP, (*R).X.as_mut_ptr());
        fpcopy(
            &Montgomery_one as *const [uint64_t; 12] as *mut digit_t as *const digit_t,
            (*R).Z.as_mut_ptr() as *mut digit_t,
        );
        fpzero((*R).Z[1 as libc::c_int as usize].as_mut_ptr());
        // Main loop
        i = 0 as libc::c_int;
        while i < nbits {
            bit = (*m.offset((i >> 6 as libc::c_int) as isize)
                >> (i & 64 as libc::c_int - 1 as libc::c_int)
                & 1 as libc::c_int as libc::c_ulong) as libc::c_int;
            swap = bit ^ prevbit;
            prevbit = bit;
            mask = (0 as libc::c_int as libc::c_ulong).wrapping_sub(swap as digit_t);
            swap_points(R, R2.as_mut_ptr(), mask);
            xDBLADD(
                R0.as_mut_ptr(),
                R2.as_mut_ptr(),
                (*R).X.as_mut_ptr() as *const felm_t,
                A24.as_mut_ptr() as *const felm_t,
            );
            fp2mul_mont(
                (*R2.as_mut_ptr()).X.as_mut_ptr() as *const felm_t,
                (*R).Z.as_mut_ptr() as *const felm_t,
                (*R2.as_mut_ptr()).X.as_mut_ptr(),
            );
            i += 1
        }
        swap = 0 as libc::c_int ^ prevbit;
        mask = (0 as libc::c_int as libc::c_ulong).wrapping_sub(swap as digit_t);
        swap_points(R, R2.as_mut_ptr(), mask);
    }
}
