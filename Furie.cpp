#include <iostream>
#include <math.h>
#include <vector>
#include <cstdlib>
using namespace std;

class Complex {
public:
    double Re, Im;
    Complex();
    Complex(double R);
    Complex(double R, double I);
    Complex operator +(const Complex &operand2) const;
    Complex operator *(double operand2) const;    
    Complex operator -(const Complex &operand2) const;
    Complex operator /(const double operand2) const;
    Complex operator ^(const double degree) const;
    friend std::ostream &operator<<(std::ostream &cout_object, const Complex operand);
};

Complex::Complex() {
    Re = 0;
    Im = 0;
}

Complex::Complex(double Re) {
    Complex::Re = Re;
    Im = 0;
}

Complex::Complex(double Re, double Im) {
    Complex::Re = Re;
    Complex::Im = Im;
}

Complex Complex::operator +(const Complex &operand2) const {
    return Complex(Re + operand2.Re, Im + operand2.Im);
}

Complex Complex::operator -(const Complex &operand2) const {
    return Complex(Re - operand2.Re, Im - operand2.Im);
}


Complex Complex::operator /(double operand2) const {
    return Complex(Re / operand2, Im / operand2);
}

Complex Complex::operator *(double operand2) const {
    return Complex(Re * operand2, Im * operand2);
}

Complex Complex::operator ^(double degree) const {
    return Complex(Re * cos(degree ) - Im * sin(degree), Im * cos(degree ) + Re * sin(degree));
}

std::ostream &operator << (std::ostream &cout_object, const Complex operand) {
    cout_object << operand.Re;
    if (operand.Im < 0) {
        cout_object << operand.Im;
    } else {
        cout_object << "+" << operand.Im;
    }
    cout_object << "i" << endl;
    return cout_object;
}

std::vector<Complex> f2s(const std::vector<Complex> &f) {
    std::vector<Complex> c;
    int n = f.size();
    c.resize(n);
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            c[k] = c[k] + (f[j] ^ (((-2) * M_PI * k * j) / n));
        }
        c[k] = c[k] / n;
    }
    return c;
}

/*std::vector<Complex> ff2s(const std::vector<Complex> &f) {
    std::vector<Complex> c, f_even, f_odd, c_even, c_odd;
    int n = f.size();
    cout << n << endl;
    int n_odd, n_even;
    if (n == 1) {
        return f;
    }        
    c.resize(n);
    if (n % 2 == 0) {
        n_odd = n / 2;
    } else {
        n_odd = n / 2 + 1;
    }
    f_even.resize(n_odd);
    n_even = n - n_odd;
    f_odd.resize( n_even);
    for (int i = 0, j = 0; j < n_even; i = i + 2, j++) {
        f_even[j] = f[i];
    }
    for (int i = 1, j = 0; j < n_odd ; i = i + 2, j++) {
        f_odd[j] = f[i];
    }
    if (n_even % 2 == 0) {
        f_even = ff2s(f_even);
    }
    else {
        f_even.push_back(0);
        
        f_even = ff2s(f_even);
        for (int i = 0; i < f_even.size(); i++) {
            cout << "EVEN" << f_even[i] << endl;
        }
    }
    if (n_odd % 2 == 0){
        f_odd = ff2s(f_odd);
    }
    else {
        f_odd.push_back(0);
        ff2s(f_odd);
        for (int i = 0; i < f_odd.size(); i++) {
            cout << "ODD" << f_odd[i] << endl;
        }
    }
    c_odd.resize(f_odd.size());
    c_even.resize(f_even.size());
    for (int i = 0; i < n_odd; i++) {
        c_odd[i] = c_odd[i] - (f_odd[i] ^ (((-2) * M_PI * i) / (n * 2)));  
        //c_odd[i] = c_odd[i] / (n);
    }
    for (int i = 0; i < n_even; i++) {
        c_even[i] = (f_even[i] ^ (((-2) * M_PI * i) / (n * 2)));
        //c_even[i] = c_even[i] / (n);
    }
    
    c_odd.insert( c_odd.end(), c_even.begin(), c_even.end() );
    return c_odd;
}*/

/*std::vector<Complex> ff2s(const std::vector<Complex> &f) {
    std::vector<Complex> c, f_even, f_odd, c_even, c_odd;
    int n = f.size();
    int n_odd, n_even;
    if (n == 1) {
        return f;
    }        
    c.resize(n);
    if (n % 2 != 0) {
        f.push_back(0);
        n++;
    }

    n_odd = n / 2;
    f_even.resize(n_odd);
    n_even = n - n_odd;
    f_odd.resize( n_even);
    for (int i = 0, j = 0; j < n_even; i++, j++) {
        f_even[j] = f[i];
    }
    for (int i = n_even, j = 0; j < n ; i = i++, j++) {
        f_odd[j] = f[i];
    }
    if (n_even % 2 == 0) {
        f_even = ff2s(f_even);
    }
    else {
        f_even.push_back(0);
        
        f_even = ff2s(f_even);
        for (int i = 0; i < f_even.size(); i++) {
            cout << "EVEN" << f_even[i] << endl;
        }
    }
    if (n_odd % 2 == 0){
        f_odd = ff2s(f_odd);
    }
    else {
        f_odd.push_back(0);
        ff2s(f_odd);
        for (int i = 0; i < f_odd.size(); i++) {
            cout << "ODD" << f_odd[i] << endl;
        }
    }
    for (int i = 0; i < n / 2 - 1; i++) {
        Complex t = x[i];
        x[i] = t + (x[i + (n / 2)] ^ ((2 * M_PI * i) / n))  
    }
    for (int i = 0)
}*/

void fft_rec(std::vector<Complex> &vec_f) {
    int n = vec_f.size(), k = n / 2;
    if (n == 1) {
        return;
    }
    std::vector<Complex> even, odd;
    even.resize(k);
    odd.resize(k);
    for (int i = 0; i < k; i++) {
        even[i] = vec_f[i << 1];
        odd[i] = vec_f[(i << 1) + 1];
    }
    fft_rec(even);
    fft_rec(odd);
    for (int i = 0; i < k; i++) {
        vec_f[i] = even[i] + (odd[i] ^ (-2 * M_PI * i / n));
        vec_f[i + k] = even[i] - (odd[i] ^ (-2 * M_PI * i / n));
    }
}

void fft(std::vector<Complex> &vec_f) {
    int size = vec_f.size(), n = 1;
    while (size > n) {
        n *= 2;
    }
    vec_f.resize(n);
    for (int i = size; i < n; i++) {
        vec_f[i] = 0;
    }
    fft_rec(vec_f);
}

void fft_rec_i(std::vector<Complex> &vec_f) {
    int n = vec_f.size(), k = n / 2;
    if (n == 1) {
        return;
    }
    std::vector<Complex> even, odd;
    even.resize(k);
    odd.resize(k);
    for (int i = 0; i < k; i++) {
        even[i] = vec_f[i * 2];
        odd[i] = vec_f[(i * 2) + 1];
    }
    fft_rec(even);
    fft_rec(odd);
    for (int i = 0; i < k; i++) {
        vec_f[i] = even[i] + (odd[i] ^ (2 * M_PI * i / n));
        vec_f[i + k] = even[i] - (odd[i] ^ (2 * M_PI * i / n));
    }
}

void fft_i(std::vector<Complex> &vec_f) {
    int size = vec_f.size(), n = 1;
    while (size > n) {
        n *= 2;
    }
    vec_f.resize(n);
    for (int i = size; i < n; i++) {
        vec_f[i] = 0;
    }
    fft_rec(vec_f);
}

std::vector <Complex> s2f(const std::vector<Complex> &c) {
    std::vector<Complex> f;
    int n = c.size();
    f.resize(n);
    for ( int j = 0; j < n; j++) {
        for(int k = 0; k < n; k++) {
            f[j] = f[j] + (c[k] ^ ((2 * M_PI * k * j) / n));
        }
    }
    return f;
}

int main() {
	Complex number;
	std::vector<Complex> vec_f = {1 , 6 , 2 , 5 , 3 , 4}, Spectr, vec_check;
    cout << "Start F:" << endl;
	for (auto &num: vec_f) {
		cout << num << endl;
	}
	Spectr = f2s(vec_f);
	vec_check = s2f(Spectr);
    cout << "SPECTR of F with DFT:" << endl;
	for (auto &num: Spectr) {
		cout << num << endl;;
	}
    cout << "F with DFT:" << endl;
	for (auto &num: vec_check) {
		cout << num << endl;
	}
	cout << "SPECTR with FFT :\n";
	fft(vec_f);
	for (auto &num: vec_f) {
		cout << num << endl;
	}
	vec_f  = s2f(vec_f);
    cout << "F with FFT" << endl;
    int n = vec_f.size();
    while(n & (n - 1) != 0) {
        n++;
    }
	for (auto &num: vec_f) {
		cout <<  num / n << endl;;
	}
	return 0;
}