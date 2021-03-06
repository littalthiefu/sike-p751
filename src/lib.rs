mod ec_isogeny;
mod fp_generic;
mod fpx;
mod sidh;
mod util;

use std::convert::From;

type uint64_t = libc::c_ulong;
type uint128_t = u128;
type digit_t = uint64_t;
type felm_t = [u64; 12];
type f2elm_t = [felm_t; 2];
type dfelm_t = [u64; 24];

const felm: felm_t = [0u64; 12];
const f2elm: f2elm_t = [felm; 2];

struct F2ELM([u8; 192]);
struct U64(Vec<u64>);
struct PointProjRaw([u64; 48]);

#[derive(Copy, Clone)]
pub struct PointProj {
    pub x: f2elm_t,
    pub z: f2elm_t,
}

impl PointProj {
    pub fn new() -> PointProj {
        PointProj {
            x: [[0u64; 12]; 2],
            z: [[0u64; 12]; 2],
        }
    }
}

type point_proj_t = [point_proj; 1];

#[derive(Copy, Clone)]
pub struct point_proj {
    X: f2elm_t,
    Z: f2elm_t,
}

impl From<PointProj> for PointProjRaw {
    fn from(a: PointProj) -> PointProjRaw {
        let mut array = [0u64; 48];

        array[..12].clone_from_slice(&a.x[0]);
        array[12..24].clone_from_slice(&a.x[1]);
        array[24..36].clone_from_slice(&a.z[0]);
        array[36..48].clone_from_slice(&a.z[1]);

        PointProjRaw(array)
    }
}

impl From<[u64; 48]> for PointProj {
    fn from(a: [u64; 48]) -> PointProj {
        let mut point = PointProj::new();

        point.x[0].clone_from_slice(&a[..12]);
        point.x[1].clone_from_slice(&a[12..24]);
        point.z[0].clone_from_slice(&a[24..36]);
        point.z[1].clone_from_slice(&a[36..48]);

        point
    }
}

impl From<f2elm_t> for F2ELM {
    fn from(a: f2elm_t) -> F2ELM {
        let mut array = [0u8; 192];
        let mut j = 0;

        for i in a.iter() {
            for i in i.iter() {
                let a = i.to_ne_bytes();
                array[j] = a[0];
                array[j + 1] = a[1];
                array[j + 2] = a[2];
                array[j + 3] = a[3];
                array[j + 4] = a[4];
                array[j + 5] = a[5];
                array[j + 6] = a[6];
                array[j + 7] = a[7];

                j += 8;
            }
        }

        F2ELM(array)
    }
}

impl From<F2ELM> for f2elm_t {
    fn from(a: F2ELM) -> f2elm_t {
        let mut array = [[0u64; 12]; 2];
        let mut i = 0;
        let mut k = 0;

        for j in (0..a.0.len()).step_by(8) {
            let mut byte = [0u8; 8];
            byte.clone_from_slice(&a.0[j..j + 8]);

            array[k][i] = u64::from_ne_bytes(byte);

            i += 1;

            if i == 12 && k != 1 {
                i = 0;
                k = 1;
            } else if i == 12 && k == 1 {
                break;
            }
        }

        array
    }
}

impl From<&[u8]> for U64 {
    fn from(a: &[u8]) -> U64 {
        let mut array = Vec::new();

        for j in (0..a.len()).step_by(8) {
            let mut byte = [0u8; 8];
            byte.clone_from_slice(&a[j..j + 8]);

            array.push(u64::from_ne_bytes(byte));
        }

        U64(array)
    }
}

const LOG2RADIX: usize = 6;
const RADIX: usize = 64;
const NWORDS_FIELD: usize = 12;

const NBITS_FIELD: usize = 751;
const MAXBITS_FIELD: usize = 768;
const MAXWORDS_FIELD: usize = (MAXBITS_FIELD + RADIX - 1) / RADIX;
const NWORDS64_FIELD: usize = (NBITS_FIELD + 63) / 64;
const NBITS_ORDER: usize = 384;
const NWORDS_ORDER: usize = (NBITS_ORDER + RADIX - 1) / RADIX;
const NWORDS64_ORDER: usize = (NBITS_ORDER + 63) / 64;
const MAXBITS_ORDER: usize = NBITS_ORDER;
const MAXWORDS_ORDER: usize = (MAXBITS_ORDER + RADIX - 1) / RADIX;
const ALICE: usize = 0;
const BOB: usize = 1;
const OALICE_BITS: usize = 372;
const OBOB_BITS: usize = 379;
const OBOB_EXPON: usize = 239;
const MASK_ALICE: usize = 0x0F;
const MASK_BOB: usize = 0x03;
const PARAM_A: usize = 6;
const PARAM_C: usize = 1;
const MAX_INT_POINTS_ALICE: usize = 8;
const MAX_INT_POINTS_BOB: usize = 10;
const MAX_Alice: usize = 186;
const MAX_Bob: usize = 239;
const MSG_BYTES: usize = 32;
const SECRETKEY_A_BYTES: usize = (OALICE_BITS + 7) / 8;
const SECRETKEY_B_BYTES: usize = (OBOB_BITS - 1 + 7) / 8;
const FP2_ENCODED_BYTES: usize = 2 * ((NBITS_FIELD + 7) / 8);

const p751: [u64; NWORDS64_FIELD] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xEEAFFFFFFFFFFFFF,
    0xE3EC968549F878A8,
    0xDA959B1A13F7CC76,
    0x084E9867D6EBE876,
    0x8562B5045CB25748,
    0x0E12909F97BADC66,
    0x00006FE5D541F71C,
];
const p751p1: [u64; NWORDS64_FIELD] = [
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0xEEB0000000000000,
    0xE3EC968549F878A8,
    0xDA959B1A13F7CC76,
    0x084E9867D6EBE876,
    0x8562B5045CB25748,
    0x0E12909F97BADC66,
    0x00006FE5D541F71C,
];
const p751x2: [u64; NWORDS64_FIELD] = [
    0xFFFFFFFFFFFFFFFE,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xDD5FFFFFFFFFFFFF,
    0xC7D92D0A93F0F151,
    0xB52B363427EF98ED,
    0x109D30CFADD7D0ED,
    0x0AC56A08B964AE90,
    0x1C25213F2F75B8CD,
    0x0000DFCBAA83EE38,
];

const Alice_order: [u64; NWORDS64_ORDER] = [
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0010000000000000,
];

const Bob_order: [u64; NWORDS64_ORDER] = [
    0xC968549F878A8EEB,
    0x59B1A13F7CC76E3E,
    0xE9867D6EBE876DA9,
    0x2B5045CB25748084,
    0x2909F97BADC66856,
    0x06FE5D541F71C0E1,
];

const A_gen: [u64; 6 * NWORDS64_FIELD] = [
    0x884F46B74000BAA8,
    0xBA52630F939DEC20,
    0xC16FB97BA714A04D,
    0x082536745B1AB3DB,
    0x1117157F446F9E82,
    0xD2F27D621A018490,
    0x6B24AB523D544BCD,
    0x9307D6AA2EA85C94,
    0xE1A096729528F20F,
    0x896446F868F3255C,
    0x2401D996B1BFF8A5,
    0x00000EF8786A5C0A,
    0xAEB78B3B96F59394,
    0xAB26681E29C90B74,
    0xE520AC30FDC4ACF1,
    0x870AAAE3A4B8111B,
    0xF875BDB738D64EFF,
    0x50109A7ECD7ED6BC,
    0x4CC64848FF0C56FB,
    0xE617CB6C519102C9,
    0x9C74B3835921E609,
    0xC91DDAE4A35A7146,
    0x7FC82A155C1B9129,
    0x0000214FA6B980B3,
    0x0F93CC38680A8CA9,
    0x762E733822E7FED7,
    0xE549F005AC0ADB67,
    0x94A71FDD2C43A4ED,
    0xD48645C2B04721C5,
    0x432DA1FE4D4CA4DC,
    0xBC99655FAA7A80E8,
    0xB2C6D502BCFD4823,
    0xEE92F40CA2EC8BDB,
    0x7B074132EFB6D16C,
    0x3340B46FA38A7633,
    0x0000215749657F6C,
    0xECFF375BF3079F4C,
    0xFBFE74B043E80EF3,
    0x17376CBE3C5C7AD1,
    0xC06327A7E29CDBF2,
    0x2111649C438BF3D4,
    0xC1F9298261BA2E97,
    0x1F9FECE869CFD1C2,
    0x01A39B4FC9346D62,
    0x147CD1D3E82A3C9F,
    0xDE84E9D249E533EE,
    0x1C48A5ADFB7C578D,
    0x000061ACA0B82E1D,
    0x1600C525D41059F1,
    0xA596899A0A1D83F7,
    0x6BFDEED6D2B23F35,
    0x5C7E707270C23910,
    0x276CA1A4E8369411,
    0xB193651A602925A0,
    0x243D239F1CA1F04A,
    0x543DC6DA457860AD,
    0xCDA590F325181DE9,
    0xD3AB7ACFDA80B395,
    0x6C97468580FDDF7B,
    0x0000352A3E5C4C77,
    0x9B794F9FD1CC3EE8,
    0xDB32E40A9B2FD23E,
    0x26192A2542E42B67,
    0xA18E94FCA045BCE7,
    0x96DC1BC38E7CDA2D,
    0x9A1D91B752487DE2,
    0xCC63763987436DA3,
    0x1316717AACCC551D,
    0xC4C368A4632AFE72,
    0x4B6EA85C9CCD5710,
    0x7A12CAD582C7BC9A,
    0x00001C7E240149BF,
];

const B_gen: [u64; 6 * NWORDS64_FIELD] = [
    0x85691AAF4015F88C,
    0x7478C5B8C36E9631,
    0x7EF2A185DE4DD6E2,
    0x943BBEE46BEB9DC7,
    0x1A3EC62798792D22,
    0x791BC4B084B31D69,
    0x03DBE6522CEA17C4,
    0x04749AA65D665D83,
    0x3D52B5C45EF450F3,
    0x0B4219848E36947D,
    0xA4CF7070466BDE27,
    0x0000334B1FA6D193,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x8E7CB3FA53211340,
    0xD67CE54F7A05EEE0,
    0xFDDC2C8BCE46FC38,
    0x08587FAE3110DF1E,
    0xD6B8246FA22B058B,
    0x4DAC3ACC905A5DBD,
    0x51D0BF2FADCED3E8,
    0xE5A2406DF6484425,
    0x907F177584F671B8,
    0x4738A2FFCCED051C,
    0x2B0067B4177E4853,
    0x00002806AC948D3D,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0xB56457016D1D6D1C,
    0x03DECCB38F39C491,
    0xDFB910AC8A559452,
    0xA9D0F17D1FF24883,
    0x8562BBAF515C248C,
    0x249B2A6DDB1CB67D,
    0x3131AF96FB46835C,
    0xE10258398480C3E1,
    0xEAB5E2B872D4FAB1,
    0xB71E63875FAEB1DF,
    0xF8384D4F13757CF6,
    0x0000361EC9B09912,
    0x58C967899ED16EF4,
    0x81998376DC622A4B,
    0x3D1C1DCFE0B12681,
    0x9347DEBB953E1730,
    0x9ABB344D3A82C2D7,
    0xE4881BD2820552B2,
    0x0037247923D90266,
    0x2E3156EDB157E5A5,
    0xF86A46A7506823F7,
    0x8FE5523A7B7F1CFC,
    0xFA3CFFA38372F67B,
    0x0000692DCE85FFBD,
];

const Montgomery_R2: [u64; NWORDS64_FIELD] = [
    0x233046449DAD4058,
    0xDB010161A696452A,
    0x5E36941472E3FD8E,
    0xF40BFE2082A2E706,
    0x4932CCA8904F8751,
    0x1F735F1F1EE7FC81,
    0xA24F4D80C1048E18,
    0xB56C383CCDB607C5,
    0x441DD47B735F9C90,
    0x5673ED2C6A6AC82A,
    0x06C905261132294B,
    0x000041AD830F1F35,
];

const Montgomery_one: [u64; NWORDS64_FIELD] = [
    0x00000000000249ad,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x8310000000000000,
    0x5527b1e4375c6c66,
    0x697797bf3f4f24d0,
    0xc89db7b2ac5c4e2e,
    0x4ca4b439d2076956,
    0x10f7926c7512c7e9,
    0x00002d5b24bce5e2,
];

const strat_Alice: [u32; MAX_Alice - 1] = [
    80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1,
    1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12,
    7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1,
    1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5,
    3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8,
    4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1,
];

const strat_Bob: [u32; MAX_Bob - 1] = [
    112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1,
    1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1,
    1, 2, 1, 1, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2,
    1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2,
    1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4,
    2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3,
    2, 1, 1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1,
    1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1,
];
