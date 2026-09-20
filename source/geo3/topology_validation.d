module geo3.topology_validation;

static import euclid_core.ring_validation;

alias RingValidationIssue =
    euclid_core.ring_validation.RingValidationIssue;

alias RingValidationResult =
    euclid_core.ring_validation.RingValidationResult;

@safe unittest
{
    RingValidationResult result;

    assert(result.valid);

    assert(
        result.issue ==
        RingValidationIssue.none
    );
}
