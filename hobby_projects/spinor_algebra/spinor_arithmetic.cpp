#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

struct spinor;
struct spinor_product;
struct spinor_sum;
struct spinor_fraction;

template <typename T, typename = void> 
struct NestingDepth { static constexpr int value = 0; };

struct spinor{
    double factor = 1;
    char bracket; 
    int i;
    int j;

    spinor(int i,int j, char bracket) {
        if (i > j){
            std::swap(i,j);
            factor *= -1;
        }
        this->i = i;
        this->j = j;
        this->bracket = bracket;
    }

    spinor_sum operator+=(const spinor s);
    spinor_sum operator-=(const spinor s);
    spinor_product operator*=(const spinor s);
    spinor_product operator/=(const spinor s);

    spinor(){}
};

struct spinor_product{
    std::vector<spinor> numerator;
    std::vector<spinor> denominator;
    double factor = 1;

    spinor_product(const spinor& s1, const spinor& s2) : numerator{s1, s2}, factor(1) {
        factor = s1.factor * s2.factor;
    }

    spinor_product(const spinor& s) : numerator{s}, factor(s.factor) {
        numerator[0].factor = 1;
    }

    spinor_product(){
        factor = 1;
    }

    spinor_sum operator+=(const spinor_product p);
    spinor_sum operator-=(const spinor_product p);
    spinor_product operator*=(const spinor_product other);
    spinor_product operator/=(const spinor_product other);
};

struct spinor_sum{
    std::vector<spinor_product> terms;
    bool is_null = false;

    spinor_sum(){
    }

    spinor_sum(const spinor_product &p1, const spinor_product &p2){
        this->terms.push_back(p1);
        this->terms.push_back(p2);
    }
        
    spinor_sum(const spinor &s1, const spinor &s2){
        spinor_product p1(s1);
        spinor_product p2(s2);

        this->terms.push_back(p1);
        this->terms.push_back(p2);
    }

    spinor_sum(const spinor_product &p){
        this->terms.push_back(p);
    }

    spinor_sum(const spinor &s){
        spinor_product p(s);
        this->terms.push_back(p);
    }

    spinor_sum& operator=(const spinor_sum& other) {
        if (this != &other) {
            this->terms = other.terms;
        }
        return *this;
    }

    spinor_sum operator+=(const spinor_sum other);
    spinor_sum operator-=(const spinor_sum other);
    spinor_sum operator*=(const spinor_sum other);
    spinor_fraction operator/=(const spinor_sum other);
};

struct spinor_fraction{
    spinor_sum upper;
    spinor_sum lower;

    spinor_fraction(){}

    spinor_fraction(const spinor &s){
        spinor_sum temp(s);
        this->upper = temp; 
    }

    spinor_fraction(const spinor_product &p){
        spinor_sum temp(p);
        this->upper = temp; 
    }

    spinor_fraction(const spinor_sum &sum){
        this->upper = sum; 
    }

    spinor_fraction(const spinor_sum &sum1, const spinor_sum &sum2){
        this->upper = sum1;
        this->lower = sum2;  
    }

    spinor_fraction& operator=(const spinor_fraction& other) {
        if (this != &other) {
            upper = other.upper;
            lower = other.lower;
        }
        return *this;
    }

    spinor_fraction operator+=(const spinor_fraction &frac);
    spinor_fraction operator-=(const spinor_fraction &frac);
    spinor_fraction operator*=(const spinor_fraction &frac);
    spinor_fraction operator/=(const spinor_fraction &frac);
};

template <typename T> 
requires requires(T t) { t.upper; } || requires(T t) { t.terms; } || requires(T t) { t.numerator; }
struct NestingDepth<T>{
    static constexpr int value = 1 + [](){
        if constexpr (requires(T t){t.upper;}){
            return NestingDepth<decltype(T::upper)>::value;
        }
        else if constexpr (requires(T t){t.terms;}){
            using InnerType = typename decltype(T::terms)::value_type;
            return NestingDepth<InnerType>::value;
        }
        else{
            using InnerType = typename decltype(T::numerator)::value_type;
            return NestingDepth<InnerType>::value;
        }
    }();
};


template<typename T, typename U>
auto operator+(T spin1, U spin2){
    constexpr int i = NestingDepth<T>::value;
    constexpr int j = NestingDepth<U>::value;

    if constexpr(i > j){
        spin1 += T(spin2);
        return spin1;
    }
    else{
        spin2 += U(spin1);
        return spin2;
    }
}

template<typename T, typename U>
auto operator-(T spin1, U spin2){
    constexpr int i = NestingDepth<T>::value;
    constexpr int j = NestingDepth<U>::value;

    if constexpr(i > j){
        spin1 -= T(spin2);
        return spin1;
    }
    else{
        spin2 -= U(spin1);
        return spin2;
    }
}

template<typename T, typename U>
auto operator*(T spin1, U spin2){
    constexpr int i = NestingDepth<T>::value;
    constexpr int j = NestingDepth<U>::value;

    if constexpr(i > j){
        spin1 *= T(spin2);
        return spin1;
    }
    else{
        spin2 *= U(spin1);
        return spin2;
    }
}

spinor_sum spinor::operator+=(const spinor s){
    return spinor_sum(*this, s);
}

spinor_sum spinor::operator-=(const spinor s){
    spinor reversed = s;
    reversed.factor *= -1;
    return spinor_sum(*this, reversed);
}

spinor_product spinor::operator*=(const spinor s){
    spinor_product result = spinor_product(*this);
    result.numerator.push_back(s);
    return result;
}

spinor_product spinor::operator/=(const spinor s){
    spinor_product result = spinor_product(*this);
    result.factor /= s.factor;
    spinor temp = s;
    temp.factor = 1;
    result.denominator.push_back(temp);
    return result;
}

spinor_sum spinor_product::operator+=(const spinor_product p){
    return spinor_sum(*this, p);
}

spinor_sum spinor_product::operator-=(const spinor_product p){
    spinor_product reversed = p;
    reversed.factor *= -1;
    return spinor_sum(*this, reversed);
}

spinor_product spinor_product::operator*=(const spinor_product other){
    numerator.insert(numerator.end(), other.numerator.begin(), other.numerator.end());
    return *this;
}

spinor_product spinor_product::operator/=(const spinor_product other){
    this->factor /= other.factor;
    numerator.insert(numerator.end(), other.denominator.begin(), other.denominator.end());
    numerator.insert(numerator.end(), other.numerator.begin(), other.numerator.end());
    return *this;
}

spinor_sum spinor_sum::operator+=(const spinor_sum other){
    terms.insert(terms.end(), other.terms.begin(), other.terms.end());
    return *this;
}

spinor_sum spinor_sum::operator-=(const spinor_sum other){
    for (spinor_product p : other.terms){
        p.factor *= -1;
        terms.push_back(p);
    }
    return *this;
}

spinor_sum spinor_sum::operator*=(const spinor_sum other){
    spinor_sum result; 
    result.terms.reserve(this->terms.size() * other.terms.size());

    for (const auto& p : this->terms) {
        for (const auto& q : other.terms) {
            result.terms.push_back(p * q);
        }
    }
    this->terms = std::move(result.terms);
    return *this;
}

spinor_fraction spinor_sum::operator/=(const spinor_sum other){
    spinor_fraction result(*this, other); 
    return result;
}

spinor_fraction spinor_fraction::operator+=(const spinor_fraction &frac){
    this->upper = this->upper * frac.lower + this->lower * frac.upper;
    this->lower = this->lower * frac.lower;
    return *this;
}

spinor_fraction spinor_fraction::operator-=(const spinor_fraction &frac){
    this->upper = this->upper * frac.lower - this->lower * frac.upper;
    this->lower = this->lower * frac.lower;
    return *this;
}

spinor_fraction spinor_fraction::operator*=(const spinor_fraction &frac){
    this->upper = this->upper * frac.upper;
    this->lower = this->lower * frac.lower;
    return *this;
}

spinor_fraction spinor_fraction::operator/=(const spinor_fraction &frac){
    this->upper = this->upper * frac.lower;
    this->lower = this->lower * frac.upper;
    return *this;
}
