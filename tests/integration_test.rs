use qrs_detector::sampling::*;
use qrs_detector::QrsDetector;

#[test]
fn test_simulated_signal() {
    let mut detector = QrsDetector::new::<216, 36>(720.sps());

    let mut detections = 0;
    let mut prev: Option<f32> = None;
    let mut prev2: Option<f32> = None;

    let data = include_str!("./data/aami3a.txt");

    let samples = data
        .split_terminator('\n')
        .map(|str| str.trim().parse::<f32>().unwrap())
        .collect::<Vec<_>>();

    for window in samples.windows(4) {
        // A slight moving average filtering
        let sum = window.iter().sum::<f32>();
        let avg = sum / 4.0;
        if let Some(p) = prev.replace(avg) {
            if let Some(p2) = prev2.replace(p) {
                if let Some(_) = detector.update((p2 - avg).abs()) {
                    detections += 1;
                }
            }
        }
    }

    assert_eq!(38, detections);
}

#[cfg(feature = "alloc")]
#[test]
fn test_simulated_signal_alloc() {
    let mut detector = QrsDetector::new_alloc(720.sps());

    let mut detections = 0;
    let mut prev: Option<f32> = None;
    let mut prev2: Option<f32> = None;

    let data = include_str!("./data/aami3a.txt");

    let samples = data
        .split_terminator('\n')
        .map(|str| str.trim().parse::<f32>().unwrap())
        .collect::<Vec<_>>();

    for window in samples.windows(4) {
        // A slight moving average filtering
        let sum = window.iter().sum::<f32>();
        let avg = sum / 4.0;
        if let Some(p) = prev.replace(avg) {
            if let Some(p2) = prev2.replace(p) {
                if let Some(_) = detector.update((p2 - avg).abs()) {
                    detections += 1;
                }
            }
        }
    }

    assert_eq!(38, detections);
}
