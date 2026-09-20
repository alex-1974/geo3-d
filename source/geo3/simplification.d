module geo3.simplification;

static import euclid_core.simplification;

alias douglasPeuckerWorkspaceSize =
    euclid_core.simplification.douglasPeuckerWorkspaceSize;

@safe unittest
{
    static assert(
        douglasPeuckerWorkspaceSize(5) ==
        3
    );
}
