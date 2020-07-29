pub fn addc(carry_in: u32, add_end_1: u64, add_end_2: u64, carry_out: &mut u32, sum_out: &mut u64) {
    let temp_reg = (add_end_1 as u128)
        .wrapping_add(add_end_2 as u128)
        .wrapping_add(carry_in as u128);

    *carry_out = (temp_reg >> 64 as i32) as u32;
    *sum_out = temp_reg as u64;
}

pub fn subc(
    borrow_in: u32,
    minu_end: u64,
    subtrah_end: u64,
    borrow_out: &mut u32,
    difference_out: &mut u64,
) {
    let temp_reg = (minu_end as u128)
        .wrapping_sub(subtrah_end as u128)
        .wrapping_sub(borrow_in as u128);

    *borrow_out = (temp_reg >> (16 as u64).wrapping_mul(8).wrapping_sub(1)) as u32;
    *difference_out = temp_reg as u64;
}

pub fn mul(multiplier: u64, multiplicand: u64, hi: &mut u64, lo: &mut u64) {
    let mut tempReg = (multiplier as u128).wrapping_mul(multiplicand as u128);
    *hi = (tempReg >> 64) as u64;
    *lo = tempReg as u64;
}

pub fn shiftl(high_in: u64, low_in: u64, shift: u32, shift_out: &mut u64) {
    *shift_out = high_in << shift ^ low_in >> (crate::RADIX as u32).wrapping_sub(shift);
}

pub fn fillbytes_u64(slice: &mut [&mut [u64]]) {
    for slice in slice.iter_mut() {
        for i in 0..slice.len() {
            let mut x = [0u8; 8];
            getrandom::getrandom(&mut x).unwrap();
            slice[i] = u64::from_be_bytes(x);
        }
    }
}
