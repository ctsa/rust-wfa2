pub mod aligner;

pub mod wfa2 {
    //! Re-export wfa2-sys bindings
    pub use wfa2_sys::*;
}

#[cfg(test)]
mod tests {
    use super::*;
    use aligner::{
        AlignmentScope, AlignmentStatus, MemoryModel, WFAlignerGapAffine, WFAlignerGapAffine2Pieces,
    };

    #[test]
    fn test_end_to_end() {
        let mut aligner =
            WFAlignerGapAffine::new(6, 4, 2, AlignmentScope::Alignment, MemoryModel::MemoryLow);

        let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusSuccessful);
        assert_eq!(aligner.score(), -24);
        assert_eq!(aligner.cigar(), b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM");
        let (a, b, c) = aligner.matching(pattern, text);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "TCTTTACTCGCGCGTT-GGAGAAATACAATAGT\n|||||||| ||||||| ||||||||||||||||\nTCTATACT-GCGCGTTTGGAGAAATAAAATAGT"
        );
    }

    #[test]
    fn test_ends_free() {
        let mut aligner = WFAlignerGapAffine::new_with_match(
            -1,
            3,
            2,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"CGCGTTTGGAGAA";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            text,
            pattern_size,
            pattern_size,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusSuccessful);
        assert_eq!(aligner.score(), 18);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(aligner.cigar(), b"IIIIIIIIIMMMMMMMMMMMMMIIIIIIIIII");

        let pattern = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let text = b"CGCGTTTGGAGAA";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            text,
            pattern_size,
            pattern_size,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusSuccessful);
        assert_eq!(aligner.score(), 18);
        assert_eq!(aligner.cigar(), b"DDDDDDDDDMMMMMMMMMMMMMDDDDDDDDDD");
    }

    #[test]
    fn test_end_to_end_affine2() {
        let mut aligner = WFAlignerGapAffine2Pieces::new_with_match(
            -1,
            3,
            3,
            3,
            10,
            0,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusSuccessful);
        assert_eq!(aligner.score(), 1);
        assert_eq!(
            std::str::from_utf8(aligner.cigar().as_slice()).unwrap(),
            "MMMMMMIIIIIIIIIIIIIIIIIIIIIMMMMM"
        );
    }
}
