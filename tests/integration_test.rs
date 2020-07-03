use if_chain::if_chain;
use qrs_detector::sampling::*;
use qrs_detector::typenum::*;
use qrs_detector::QrsDetector;

#[test]
fn test_simulated_signal() {
    let mut detector: QrsDetector<U216, U36> = QrsDetector::new(720.sps());

    let mut detections = 0;
    let mut prev: Option<f32> = None;
    let mut prev2: Option<f32> = None;

    let data = include_str!("./data/aami3a.txt");

    for sample_str in data.split_terminator('\n') {
        let sample: f32 = sample_str.trim().parse().unwrap();
        if_chain! {
            if let Some(p) = prev.replace(sample);
            if let Some(p2) = prev2.replace(p);
            if let Some(_) = detector.update((p2 - sample).abs());
            then {
                detections += 1;
            }
        }
    }

    assert_eq!(50, detections);
}
