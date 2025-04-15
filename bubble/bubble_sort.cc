/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void bubble_sort(int A[], int size) {
    for (int i = 1; i <= size - 1; i++) {
        for (int j = size; j >= i + 1; --j) {
            if (A[j] < A[j - 1]) {
                int temp = A[j - 1];
                A[j - 1] = A[j];
                A[j] = temp;
            }
        }
    }
}