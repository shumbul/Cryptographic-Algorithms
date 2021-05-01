// C++ program Miller-Rabin primality test
#include <bits/stdc++.h>
using namespace std;

// Utility function to do modular exponentiation.
// It returns (x^y) % p
int power(int x, unsigned int y, int p)
{
    int res = 1;      // Initialize result
    x = x % p;  // Update x if it is more than or
    // equal to p
    while (y > 0)
    {
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res * x) % p;

        // y must be even now
        y = y >> 1; // y = y/2
        x = (x * x) % p;
    }
    return res;
}

// This function is called for all k trials. It returns
// false if n is composite and returns false if n is
// probably prime.
// d is an odd number such that  d*2<sup>r</sup> = n-1
// for some r >= 1
bool miillerTest(int d, int n)
{
    cout << " ### Miller Rabin test for d = " << d << ", n = " << n << " ###\n";
    cout << " STEP 1: \nPick a random number in [2..n-2]\n Corner cases make sure that n > 4\n";
    int a = 2 + rand() % (n - 4);

    cout << " Chosen random number is " << a << "\n";
    cout << "\nSTEP 2: \nCompute a^d % n \n";
    int x = power(a, d, n);
    cout << " = " << x << "\n\n";

    if (x == 1  || x == n - 1)
        return true;

    // Keep squaring x while one of the following doesn't
    // happen
    // (i)   d does not reach n-1
    // (ii)  (x^2) % n is not 1
    // (iii) (x^2) % n is not n-1
    while (d != n - 1)
    {
        cout << "At d = " << d << "\tx = " << x << "\n";
        x = (x * x) % n;
        d *= 2;

        if (x == 1)      return false;
        if (x == n - 1)    return true;
    }

    // Return composite
    return false;
}

// It returns false if n is composite and returns true if n
// is probably prime.  k is an input parameter that determines
// accuracy level. Higher value of k indicates more accuracy.
bool isPrime(int n, int k)
{
    // Corner cases
    if (n <= 1 || n == 4)  return false;
    if (n <= 3) return true;

    // Find r such that n = 2^d * r + 1 for some r >= 1
    int d = n - 1;
    while (d % 2 == 0)
        d /= 2;

    // Iterate given nber of 'k' times
    for (int i = 0; i < k; i++) {
        cout << "*** TEST ***" << i + 1 << " :\n\n";
        if (!miillerTest(d, n))
            return false;
        cout << " test " << i + 1 << " PASSED!\n";
    }

    return true;
}

// Driver program
int main()
{
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif
    int n;
    cin >> n;
    int k = 4;
    if (isPrime(n, k))
        cout << n << " is prime\n";
    else
        cout << n << " is NOT prime\n";

    return 0;
}