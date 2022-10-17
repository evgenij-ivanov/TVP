#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Class.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

int main()
{
    try
    {
        freopen("input1.txt", "r", stdin);
        int n;
        cin >> n;
        vector<vector<double>> a(n, vector<double>(n));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cin >> a[i][j];
            }
        }
        Matrix<double> matrix(a);
        auto start = high_resolution_clock::now();
        auto inversedMatrix = matrix.inverse();
        auto stop = high_resolution_clock::now();
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << inversedMatrix->at(i, j) << ' ';
            }
            cout << endl;
        }
        delete inversedMatrix;
        auto duration = duration_cast<nanoseconds>(stop - start);
        cout << duration.count() << endl;
        start = high_resolution_clock::now();
        auto parallelInversedMatrix = matrix.parallelInverse();
        stop = high_resolution_clock::now();
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << parallelInversedMatrix->at(i, j) << ' ';
            }
            cout << endl;
        }
        delete parallelInversedMatrix;
        duration = duration_cast<nanoseconds>(stop - start);
        cout << duration.count() << endl;
    }
    catch (const std::invalid_argument& e)
    {
        cerr << e.what() << endl;
    }
}
