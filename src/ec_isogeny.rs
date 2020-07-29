use crate::{f2elm_t, fpx, PointProj};

pub fn xDBL(p: &PointProj, q: &mut PointProj, a24plus: &f2elm_t, c24: &f2elm_t) {
    let mut t0: f2elm_t = [[0u64; 12]; 2];
    let mut t1: f2elm_t = [[0u64; 12]; 2];

    fpx::fp2sub(&p.x, &p.z, &mut t0);
    fpx::fp2add(&p.x, &p.z, &mut t1);
    let temp = t0;
    fpx::fp2sqr_mont(&temp, &mut t0);
    let temp = t1;
    fpx::fp2sqr_mont(&temp, &mut t1);
    fpx::fp2mul_mont(c24, &t0, &mut q.z);
    fpx::fp2mul_mont(&t1, &q.z, &mut q.x);
    let temp = t1;
    fpx::fp2sub(&temp, &t0, &mut t1);
    fpx::fp2mul_mont(a24plus, &t1, &mut t0);
    let temp = q.z;
    fpx::fp2add(&temp, &t0, &mut q.z);
    let temp = q.z;
    fpx::fp2mul_mont(&temp, &t1, &mut q.z);
}

// https://discord.com/channels/244230771232079873/287491279888318465/737645131372691467
pub fn xDBLe(p: &PointProj, q: &mut PointProj, a24plus: &f2elm_t, c24: &f2elm_t, e: usize) {
    fpx::copy_words(&p.x[0], &mut q.x[0], 2 * 2 * crate::NWORDS_FIELD);

    for _ in 0..e {
        let temp = *q;
        xDBL(&temp, q, a24plus, c24);
    }
}

pub fn get_4_isog(p: &PointProj, a24plus: &mut f2elm_t, c24: &mut f2elm_t, coeff: &mut [f2elm_t]) {
    fpx::fp2sub(&p.x, &p.z, &mut coeff[1]);
    fpx::fp2add(&p.x, &p.z, &mut coeff[2]);
    fpx::fp2sqr_mont(&p.z, &mut coeff[0]);
    let temp = coeff[0];
    fpx::fp2add(&temp, &temp, &mut coeff[0]);
    fpx::fp2sqr_mont(&coeff[0], c24);
    let temp = coeff[0];
    fpx::fp2add(&temp, &temp, &mut coeff[0]);
    fpx::fp2sqr_mont(&p.x, a24plus);
    let temp = *a24plus;
    fpx::fp2add(&temp, &temp, a24plus);
    let temp = *a24plus;
    fpx::fp2sqr_mont(&temp, a24plus);
}
