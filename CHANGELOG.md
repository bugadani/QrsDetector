Unreleased
==========

 * **breaking:** Replaced `typenum` with const generics.
 * `SamplingFrequency` internal field is no longer accessible.
 * `SamplingFrequencyExt` is now implemented for `f64`.
 * Added `Thresholds`, `QrsDetector::thresholds()`.
 * Added `alloc` feature to automatically allocate F buffers.

0.2.0
==========

 * **breaking:** Remove preprocessing. Users must derive their differentiated complex leads.
 * **breaking:** Internal buffer size to calculate `F` changed from 350ms to 300ms.

0.1.2
=====

 * Small internal improvements

0.1.1
=====
 * Dependency update

0.1.0
=====
 * Initial release