#include <cmath>
#include <climits>

/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */
/* 1-indexed algorithm! */

using namespace std;

struct Node {
    int data;
    int child_index;
};


void tournament_sort(int A[], int size) {
    int leaf_count = (int)pow(2, ceil(log2(size)));
    int total_nodes = 2 * leaf_count;
    Node* tree = new Node[total_nodes];

    for (int i = 0; i < total_nodes; ++i) {
        tree[i].data = INT_MAX;
        tree[i].data = INT_MAX;
        tree[i].child_index = -1;
    }
    for (int i = 0; i < size; ++i) {
        tree[leaf_count + i].data = A[i + 1];
    }

    // Build tree from bottom-up
    for (int i = leaf_count - 1; i > 0; --i) {
        int left = 2 * i;
        int right = 2 * i + 1;
        if (tree[left].data < tree[right].data) {
            tree[i].data = tree[left].data;
            tree[i].child_index = left;
        } else {
            tree[i].data = tree[right].data;
            tree[i].child_index = right;
        }
    }

    for (int i = 1; i <= size; ++i) {
        A[i] = tree[1].data;

        int idx = 1;
        while (tree[idx].child_index != -1) {
            idx = tree[idx].child_index;
        }
        tree[idx].data = INT_MAX;
        tree[idx].data = INT_MAX;

        idx /= 2; // up to parent
        while (idx > 0) {
            int left = 2 * idx;
            int right = 2 * idx + 1;
            if (tree[left].data < tree[right].data) {
                tree[idx].data = tree[left].data;
                tree[idx].child_index = left;
            } else {
                tree[idx].data = tree[right].data;
                tree[idx].child_index = right;
            }
            idx /= 2;
        }
    }

    delete[] tree;
}