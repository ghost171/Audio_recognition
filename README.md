# Audio_recognition
## Discrete and fast Fourier transform
### Explaination
In mathematics, a Fourier transform (FT) is a mathematical transform which decomposes a function (often a function of time, or a signal) into its constituent frequencies, such as the expression of a musical chord in terms of the volumes and frequencies of its constituent notes. The term Fourier transform refers to both the frequency domain representation and the mathematical operation that associates the frequency domain representation to a function of time.

The Fourier transform of a function of time is a complex-valued function of frequency, whose magnitude (absolute value) represents the amount of that frequency present in the original function, and whose argument is the phase offset of the basic sinusoid in that frequency. The Fourier transform is not limited to functions of time, but the domain of the original function is commonly referred to as the time domain. 

There is also an inverse Fourier transform that mathematically synthesizes the original function from its frequency domain representation, as proven by the Fourier inversion theorem. 

Linear operations performed in one domain (time or frequency) have corresponding operations in the other domain, which are sometimes easier to perform. The operation of differentiation in the time domain corresponds to multiplication by the frequency, so some differential equations are easier to analyze in the frequency domain. Also, convolution in the time domain corresponds to ordinary multiplication in the frequency domain (see Convolution theorem). After performing the desired operations, transformation of the result can be made back to the time domain. Harmonic analysis is the systematic study of the relationship between the frequency and time domains, including the kinds of functions or operations that are "simpler" in one or the other, and has deep connections to many areas of modern mathematics.

Functions that are localized in the time domain have Fourier transforms that are spread out across the frequency domain and vice versa, a phenomenon known as the uncertainty principle. The critical case for this principle is the Gaussian function, of substantial importance in probability theory and statistics as well as in the study of physical phenomena exhibiting normal distribution (e.g., diffusion). The Fourier transform of a Gaussian function is another Gaussian function. Joseph Fourier introduced the transform in his study of heat transfer, where Gaussian functions appear as solutions of the heat equation.

The Fourier transform can be formally defined as an improper Riemann integral, making it an integral transform, although this definition is not suitable for many applications requiring a more sophisticated integration theory. For example, many relatively simple applications use the Dirac delta function, which can be treated formally as if it were a function, but the justification requires a mathematically more sophisticated viewpoint. The Fourier transform can also be generalized to functions of several variables on Euclidean space, sending a function of 3-dimensional 'position space' to a function of 3-dimensional momentum (or a function of space and time to a function of 4-momentum). This idea makes the spatial Fourier transform very natural in the study of waves, as well as in quantum mechanics, where it is important to be able to represent wave solutions as functions of either position or momentum and sometimes both. In general, functions to which Fourier methods are applicable are complex-valued, and possibly vector-valued. Still further generalization is possible to functions on groups, which, besides the original Fourier transform on ℝ or ℝn (viewed as groups under addition), notably includes the discrete-time Fourier transform (DTFT, group = ℤ), the discrete Fourier transform (DFT, group = ℤ mod N) and the Fourier series or circular Fourier transform (group = S1, the unit circle ≈ closed finite interval with endpoints identified). The latter is routinely employed to handle periodic functions. 
The fast Fourier transform (FFT) is an algorithm for computing the DFT. 
### Formule for DFT
f(t) = c0·1 + c1e2πti + c2e4πi + c2e6πi + ...
## FFT explanation
A fast Fourier transform (FFT) is an algorithm that computes the discrete Fourier transform (DFT) of a sequence, or its inverse (IDFT). Fourier analysis converts a signal from its original domain (often time or space) to a representation in the frequency domain and vice versa. The DFT is obtained by decomposing a sequence of values into components of different frequencies. This operation is useful in many fields, but computing it directly from the definition is often too slow to be practical. An FFT rapidly computes such transformations by factorizing the DFT matrix into a product of sparse (mostly zero) factors. As a result, it manages to reduce the complexity of computing the DFT from O(N2), which arises if one simply applies the definition of DFT, to O ( NlogN ), where N  is the data size. The difference in speed can be enormous, especially for long data sets where N may be in the thousands or millions. In the presence of round-off error, many FFT algorithms are much more accurate than evaluating the DFT definition directly or indirectly. There are many different FFT algorithms based on a wide range of published theories, from simple complex-number arithmetic to group theory and number theory.

Fast Fourier transforms are widely used for applications in engineering, music, science, and mathematics. The basic ideas were popularized in 1965, but some algorithms had been derived as early as 1805. In 1994, Gilbert Strang described the FFT as "the most important numerical algorithm of our lifetime", and it was included in Top 10 Algorithms of 20th Century by the IEEE magazine Computing in Science & Engineering.

The best-known FFT algorithms depend upon the factorization of N, but there are FFTs with O(NlogN) complexity for all N, even for prime N. Many FFT algorithms only depend on the fact that e^(−2πi/N)  is an N-th primitive root of unity, and thus can be applied to analogous transforms over any finite field, such as number-theoretic transforms. Since the inverse DFT is the same as the DFT, but with the opposite sign in the exponent and a 1/N factor, any FFT algorithm can easily be adapted for it. 
Cooley–Tukey algorithm
Main article: Cooley–Tukey FFT algorithm

By far the most commonly used FFT is the Cooley–Tukey algorithm. This is a divide and conquer algorithm that recursively breaks down a DFT of any composite size N = N1N2 into many smaller DFTs of sizes N1 and N2, along with O(N) multiplications by complex roots of unity traditionally called twiddle factors (after Gentleman and Sande, 1966).

This method (and the general idea of an FFT) was popularized by a publication of Cooley and Tukey in 1965, but it was later discovered that those two authors had independently re-invented an algorithm known to Carl Friedrich Gauss around 1805 (and subsequently rediscovered several times in limited forms).

The best known use of the Cooley–Tukey algorithm is to divide the transform into two pieces of size N/2 at each step, and is therefore limited to power-of-two sizes, but any factorization can be used in general (as was known to both Gauss and Cooley/Tukey). These are called the radix-2 and mixed-radix cases, respectively (and other variants such as the split-radix FFT have their own names as well). Although the basic idea is recursive, most traditional implementations rearrange the algorithm to avoid explicit recursion. Also, because the Cooley–Tukey algorithm breaks the DFT into smaller DFTs, it can be combined arbitrarily with any other algorithm for the DFT, such as those described below. 
![image](https://wikimedia.org/api/rest_v1/media/math/render/svg/350d2578274aae20d7db9bcccb44f1c49c8a889d)
## split_by_vad.cpp
#### What does this program?
This program supposed to dividing audio record ro few record with human voice.
It used to dividing audio record with five digits in it to 5 records and put it to the appropriate in splitted/ directories.
#### functions for vad.cpp
    def print_with_timeline(data, single_duration, units_name, row_limit):
    def get_segment_energy(data, start, end):
    def get_segments_energy(data, segment_duration):
    def get_vad_mask(data, threshold):
    def sec2samples(seconds, sample_rate):
    class Segment:
    def print_segments(segments, single_duration, units_name):
    def mask_compress(data):
#### functions and methods for Furie    
    Complex::Complex() {}
    Complex::Complex(double Re) {}
    Complex::Complex(double Re, double Im) {}
    Complex Complex::operator +(const Complex &operand2) const {}
    Complex Complex::operator -(const Complex &operand2) const {}
    Complex Complex::operator /(double operand2) const {}
    Complex Complex::operator *(double operand2) const {}
    Complex Complex::operator ^(double degree) const {}
    std::ostream &operator << (std::ostream &cout_object, const Complex operand) {}
    std::vector<Complex> f2s(const std::vector<Complex> &f) {}
    void fft_rec(std::vector<Complex> &vec_f) {}
    void fft(std::vector<Complex> &vec_f) {}
    void fft_rec_i(std::vector<Complex> &vec_f) {}
    void fft_i(std::vector<Complex> &vec_f) {}
    std::vector <Complex> s2f(const std::vector<Complex> &c) {}
