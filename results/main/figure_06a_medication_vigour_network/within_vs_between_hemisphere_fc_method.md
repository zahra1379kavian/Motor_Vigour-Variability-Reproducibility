# Within-Hemisphere vs Between-Hemisphere FC Method

This analysis uses the lateralized ROI network matrices produced by the medication-effects pipeline. ROIs must have names ending in `_L` or `_R`. For each subject/session, upper-triangle ROI edges were classified as within-hemisphere when both ROIs had the same hemisphere suffix, and between-hemisphere when one ROI was left-lateralized and the other was right-lateralized. The session-level within- and between-hemisphere values are unweighted means of the mutual-information edge weights in each class.

Medication effects were evaluated within complete subjects as ON minus OFF separately for within- and between-hemisphere edges. The primary comparison was the paired subject-level contrast (ON - OFF within-hemisphere FC) - (ON - OFF between-hemisphere FC). Positive values mean medication increased within-hemisphere connectivity more than between-hemisphere connectivity.
