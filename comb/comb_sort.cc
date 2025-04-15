/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void swap1(int A[], int p, int q) {
    int temp = A[p];
    A[p] = A[q];
    A[q] = temp;
}

void comb_sort(int A[], int size) {
    int gap = size - 2;
    const double shrink = 1.3;
    bool swapped = true;

    while (gap > 1 || swapped) {
        gap = (int)(gap / shrink);
        if (gap < 1) gap = 1;

        swapped = false;
        for (int i = 1; i + gap <= size; ++i) {
            if (A[i] > A[i + gap]) {
                swap1(A, i, i + gap);
                swapped = true;
            }
        }
    }
}
