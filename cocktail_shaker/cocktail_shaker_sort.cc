/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

void swap(int A[], int p, int q) {
    int temp = A[p];
    A[p] = A[q];
    A[q] = temp;
}

void cocktail_shaker_sort(int A[], int size) {
    int l = 1;
    int r = size;
    int it = 1;
    bool swapped = true;

    while (l != r && swapped) {
        swapped = false;
        if (it % 2 == 1) {
            for (int i = l; i < r; i++) {
                if (A[i] > A[i + 1]) {
                    swap(A, i, i+1);
                    swapped = true;
                }
            } --r;
        }
        else {
            for (int i = r; i > l; --i) {
                if (A[i] < A[i - 1]) {
                    swap(A, i, i-1);
                    swapped = true;
                }
            } l++;
        }
        it++;
    }
}