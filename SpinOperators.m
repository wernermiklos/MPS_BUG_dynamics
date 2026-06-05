function [Sx, Sy, Sz] = SpinOperators(S)
%SPINOPERATORS Spin angular-momentum operators for arbitrary spin S.
%   [Sx, Sy, Sz] = spinOperators(S) returns the (2S+1)-by-(2S+1) spin
%   operator matrices for the spin quantum number S = 0, 1/2, 1, 3/2, ...
%
%   Operators are expressed in units of hbar (hbar = 1).
%
%   Basis ordering: descending magnetic quantum number m,
%       m = S, S-1, ..., -S   ==>   Sz = diag(S, S-1, ..., -S).
%
%   Definitions used:
%       S+ |S,m> = sqrt(S(S+1) - m(m+1)) |S,m+1>,   S- = (S+)'
%       Sx = (S+ + S-)/2,   Sy = (S+ - S-)/(2i)
%
%   Example:
%       [Sx, Sy, Sz] = spinOperators(1);   % spin-1, 3x3 matrices

    % --- input validation -------------------------------------------------
    if ~isscalar(S) || S < 0 || abs(2*S - round(2*S)) > 1e-12
        error('spinOperators:badInput', ...
              'S must be a non-negative integer or half-integer (0, 1/2, 1, ...).');
    end

    N = round(2*S + 1);          % matrix dimension
    m = (S:-1:-S).';             % column of m values: S, S-1, ..., -S

    % Sz is diagonal in this basis
    Sz = diag(m);

    % Raising operator S+ : nonzero on the first superdiagonal.
    % The (j, j+1) entry uses m of column j+1, i.e. m(2:end).
    offdiag = sqrt(S*(S+1) - m(2:end).*(m(2:end) + 1));
    Sp = diag(offdiag, 1);       % S_plus
    Sm = Sp';                    % S_minus = (S_plus)^dagger

    % Cartesian components
    Sx = (Sp + Sm) / 2;
    Sy = (Sp - Sm) / (2i);
end