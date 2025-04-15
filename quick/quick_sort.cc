/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

int partition(int A[], int p, int r) {
    int x = A[r - 1]; // pivot
    int i = p - 1;
    for (int j = p; j < r - 1; ++j) {
        if (A[j] <= x) {
            i = i + 1;
            int temp = A[j];
            A[j] = A[i];
            A[i] = temp;
        }
    }
    int temp = A[i + 1]; // exchange with pivot
    A[i + 1] = A[r - 1];
    A[r - 1] = temp;
    return i + 1;
}

void quick_sort1(int A[], int p, int r) {
    if (p < r - 1) {
        int q = partition(A, p, r);
        quick_sort1(A, p, q);
        quick_sort1(A, q, r);
    }
}

void quick_sort(int A[], int size) {
    int p = 1;
    int r = size + 1;
    quick_sort1(A, p, r);
}
