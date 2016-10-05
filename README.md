


This version of dmnautilus inverts the logic in CIAO 4.6.

In the original code, the SNR is used as an upper limit.  The cells
are split until the SNR falls below the SNR limit.

In this version, the code uses the SNR as a lower limit.  The cells
are only split if each split cell has a SNR above the limit.
