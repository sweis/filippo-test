//! Drive the verifier against the Wycheproof ML-DSA test vectors.

use mldsa_verify::{params::ParamSet, verify, VerifyError, ML_DSA_44, ML_DSA_65, ML_DSA_87};
use serde::Deserialize;

#[derive(Deserialize)]
struct TestFile {
    #[serde(rename = "testGroups")]
    test_groups: Vec<TestGroup>,
}

#[derive(Deserialize)]
struct TestGroup {
    #[serde(rename = "publicKey")]
    public_key: String,
    tests: Vec<TestCase>,
}

#[derive(Deserialize)]
struct TestCase {
    #[serde(rename = "tcId")]
    tc_id: u32,
    #[serde(default)]
    comment: String,
    msg: String,
    #[serde(default)]
    ctx: Option<String>,
    sig: String,
    result: String,
    #[serde(default)]
    flags: Vec<String>,
}

fn run_file(path: &str, ps: &ParamSet) {
    let data = std::fs::read_to_string(path).expect("read test vectors");
    let file: TestFile = serde_json::from_str(&data).expect("parse test vectors");

    let mut passed = 0usize;
    let mut total = 0usize;
    let mut failures = Vec::new();

    for group in &file.test_groups {
        let pk = hex::decode(&group.public_key).expect("hex pk");
        for tc in &group.tests {
            total += 1;
            let msg = hex::decode(&tc.msg).expect("hex msg");
            let sig = hex::decode(&tc.sig).expect("hex sig");
            let ctx_bytes = tc.ctx.as_deref().map(|s| hex::decode(s).expect("hex ctx"));
            let ctx: &[u8] = ctx_bytes.as_deref().unwrap_or(&[]);

            let got = verify(ps, &pk, &msg, &sig, ctx);
            let expected_valid = tc.result == "valid";
            let ok = expected_valid == got.is_ok();
            if ok {
                passed += 1;
            } else {
                failures.push(format!(
                    "tcId {}: expected {:?}, got {:?} (flags: {:?}, comment: {})",
                    tc.tc_id,
                    tc.result,
                    got.err().map_or("Ok".to_string(), |e| format!("{:?}", e)),
                    tc.flags,
                    tc.comment
                ));
            }

            // Extra strictness: specific error classes for specific flags.
            // This sharpens mutation testing by tying specific VerifyError
            // values to the flagged failure mode.
            if tc.flags.iter().any(|f| f == "IncorrectPublicKeyLength") {
                assert_eq!(got, Err(VerifyError::BadPublicKeyLength), "tcId {}", tc.tc_id);
            }
            if tc.flags.iter().any(|f| f == "IncorrectSignatureLength") {
                assert_eq!(got, Err(VerifyError::BadSignatureLength), "tcId {}", tc.tc_id);
            }
            if tc.flags.iter().any(|f| f == "InvalidContext") {
                assert_eq!(got, Err(VerifyError::ContextTooLong), "tcId {}", tc.tc_id);
            }
            if tc.flags.iter().any(|f| f == "InvalidHintsEncoding") {
                assert_eq!(got, Err(VerifyError::MalformedHint), "tcId {}", tc.tc_id);
            }
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("{f}");
        }
    }
    assert_eq!(
        passed, total,
        "{} of {} Wycheproof tests passed for {path}",
        passed, total
    );
}

#[test]
fn wycheproof_mldsa_44() {
    run_file("testvectors/mldsa_44_verify_test.json", &ML_DSA_44);
}

#[test]
fn wycheproof_mldsa_65() {
    run_file("testvectors/mldsa_65_verify_test.json", &ML_DSA_65);
}

#[test]
fn wycheproof_mldsa_87() {
    run_file("testvectors/mldsa_87_verify_test.json", &ML_DSA_87);
}
