/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void merge(int A[], int p, int q, int r) {
    int n1 = q - p;
    int n2 = r - q;

    int L[n1 + 1];
    int R[n2 + 1];

    for (int i = 0; i < n1; ++i)
        L[i] = A[p + i];

    for (int j = 0; j < n2; ++j)
        R[j] = A[q + j];

    L[n1] = __INT32_MAX__;
    R[n2] = __INT32_MAX__;

    int i = 0, j = 0;

    for (int k = p; k < r; ++k) {
        if (L[i] <= R[j]) {
            A[k] = L[i];
            ++i;
        } else {
            A[k] = R[j];
            ++j;
        }
    }
}

void merge_sort1(int A[], int p, int r) {
    if (p < r - 1) {
        int q = (p + r) / 2;
        merge_sort1(A, p, q);
        merge_sort1(A, q, r);
        merge(A, p, q, r);
    }
}

void merge_sort(int A[], int size) {
    int p = 1;
    int r = size + 1;
    merge_sort1(A, p, r);
}
