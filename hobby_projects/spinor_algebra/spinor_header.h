#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

namespace Spinor{
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

        spinor(){};
    };

    struct spinor_product{
        std::vector<spinor> numerator;
        std::vector<spinor> denominator;
        double factor = 1;

        spinor_product(const spinor& s1, const spinor& s2) : numerator{s1, s2}, factor(1) {
        factor = s1.factor * s2.factor;
        }

        spinor_product(const spinor& s) : numerator{s}, factor(1) {
        factor = s.factor;
        }

        spinor_product(){
            factor = 1;
        }
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
        
        spinor_sum(const spinor &s1,const spinor &s2){
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
    };

    struct spinor_fraction{
        spinor_sum numerator;
        spinor_sum denominator;

        spinor_fraction(){}
        spinor_fraction(const spinor &s){
            spinor_sum temp(s);
            this->numerator = temp; 
        }

        spinor_fraction(const spinor_product &p){
            spinor_sum temp(p);
            this->numerator = temp; 
        }

        spinor_fraction(const spinor_sum &sum){
            this->numerator = sum; 
        }

        spinor_fraction(const spinor_sum &sum1, const spinor_sum &sum2){
            this->numerator = sum1;
            this->denominator = sum2;  
        }

        spinor_fraction& operator=(const spinor_fraction& other) {
            if (this != &other) {
                numerator = other.numerator;
                denominator = other.denominator;
            }
        return *this;
        }
    };

    struct factor_struct{
        spinor factor;
        int factor_count;

        factor_struct(){}

        factor_struct(const spinor &s, const int &num){
            this->factor = s;
            this->factor_count = num;
        }
    };

    struct factorized_spinor {
        std::vector<spinor_sum> factor;     
        spinor_sum remainder;  
        int original_sum_length;
        int order;

        factorized_spinor(){}
        factorized_spinor(const spinor_sum &non_factored) : remainder(non_factored) {
            this->original_sum_length = remainder.terms.size();
            this->order = original_sum_length;
        }
    };

    spinor_sum operator+(const spinor &s1, const spinor &s2);
    spinor_sum operator+(const spinor_product &p1, const spinor &s2);
    spinor_sum operator+(const spinor &s1,const spinor_product &p2);
    spinor_sum operator+(const spinor_product &p1, const spinor_product &p2);
    spinor_sum operator+(const spinor_sum &sum1, const spinor_product &p2);
    spinor_sum operator+(const spinor_product &p1, const spinor_sum &sum2);
    spinor_sum operator+(const spinor_sum &sum1, const spinor &s2);
    spinor_sum operator+(const spinor &s1,const spinor_sum &sum2);
    spinor_sum operator+(const spinor_sum &sum1, const spinor_sum &sum2);
    spinor_fraction operator+(const spinor_fraction &f1, const spinor &s2);
    spinor_fraction operator+(const spinor &s1, const spinor_fraction &f2);
    spinor_fraction operator+(const spinor_fraction &f1, const spinor_product &p2);
    spinor_fraction operator+(const spinor_product &p1, const spinor_fraction &f2);
    spinor_fraction operator+(const spinor_fraction &f1, const spinor_sum &sum2);
    spinor_fraction operator+(const spinor_sum &sum1, const spinor_fraction &f2);
    spinor_fraction operator+(const spinor_fraction &f1, const spinor_fraction &f2);

    spinor_sum operator-(const spinor &s1, const spinor &s2);
    spinor_sum operator-(const spinor_product &p1, const spinor &s2);
    spinor_sum operator-(const spinor &s1,const spinor_product &p2);
    spinor_sum operator-(const spinor_product &p1, const spinor_product &p2);
    spinor_sum operator-(const spinor_sum &sum1, const spinor_product &p2);
    spinor_sum operator-(const spinor_product &p1, const spinor_sum &sum2);
    spinor_sum operator-(const spinor_sum &sum1, const spinor &s2);
    spinor_sum operator-(const spinor &s1,const spinor_sum &sum2);
    spinor_sum operator-(const spinor_sum &sum1, const spinor_sum &sum2);
    spinor_fraction operator-(const spinor_fraction &f1, const spinor &s2);
    spinor_fraction operator-(const spinor &s1, const spinor_fraction &f2);
    spinor_fraction operator-(const spinor_fraction &f1, const spinor_product &p2);
    spinor_fraction operator-(const spinor_product &p1, const spinor_fraction &f2);
    spinor_fraction operator-(const spinor_fraction &f1, const spinor_sum &sum2);
    spinor_fraction operator-(const spinor_sum &sum1, const spinor_fraction &f2);
    spinor_fraction operator-(const spinor_fraction &f1, const spinor_fraction &f2);

    spinor_product operator*(const spinor &s1, const spinor &s2);
    spinor_product operator*(const spinor_product &p1, const spinor &s2);
    spinor_product operator*(const spinor &s1, const spinor_product &p2);
    spinor_product operator*(const spinor_product &p1, const spinor_product &p2);
    spinor operator*(const double &c, const spinor &s);
    spinor_product operator*(const double &c, const spinor_product &p);
    spinor operator*(const spinor &s,const double &c);
    spinor_product operator*(const spinor_product &p,const double &c);
    spinor_sum operator*(const spinor_product &p1,const spinor_sum &sum2);
    spinor_sum operator*(const spinor_sum &sum1,const spinor &s2);
    spinor_sum operator*(const spinor_sum &sum1,const spinor_product &p2);
    spinor_sum operator*(const spinor_sum &sum1,const spinor_sum &sum2);
    spinor_sum operator*(const spinor &s1,const spinor_sum &sum2);
    spinor_sum operator*(const spinor &sum1,const spinor_sum &s2);
    spinor_sum operator*(const double &c, const spinor_sum &sum);
    spinor_sum operator*(const spinor_sum &sum,const double &c);
    spinor_fraction operator*(const spinor_fraction &f1, const spinor &s2);
    spinor_fraction operator*(const spinor &s1, const spinor_fraction &f2);
    spinor_fraction operator*(const spinor_fraction &f1, const spinor_product &p2);
    spinor_fraction operator*(const spinor_product &p1, const spinor_fraction &f2);
    spinor_fraction operator*(const spinor_fraction &f1, const spinor_sum &sum2);
    spinor_fraction operator*(const spinor_sum &sum1, const spinor_fraction &f2);
    spinor_fraction operator*(const spinor_fraction &f1, const spinor_fraction &f2);
    spinor_fraction operator*(const double &c, const spinor_fraction &f);
    spinor_fraction operator*(const spinor_fraction &f, const double &c);

    spinor_product operator/(const spinor &s1, const spinor &s2);
    spinor_product operator/(const spinor_product &p1, const spinor &s2);
    spinor_product operator/(const spinor &s1, const spinor_product &p2);
    spinor_sum operator/(const spinor_sum &sum1,const spinor &s2);
    spinor_product operator/(const spinor_product &p1, const spinor_product &p2);
    spinor_product operator/(const double &c, const spinor &s);
    spinor_product operator/(const double &c, const spinor_product &p);
    spinor operator/(const spinor &s,const double &c);
    spinor_product operator/(const spinor_product &p,const double &c);
    spinor_fraction operator/(const spinor_product &p1,const spinor_sum &sum2);
    spinor_sum operator/(const spinor_sum &sum1,const spinor_product &p2);
    spinor_fraction operator/(const spinor_sum &sum1,const spinor_sum &sum2);
    spinor_fraction operator/(const spinor &s1,const spinor_sum &sum2);
    spinor_fraction operator/(const double &c, const spinor_sum &sum);
    spinor_sum operator/(const spinor_sum &sum,const double &c);
    spinor_fraction operator/(const spinor_fraction &f1, const spinor &s2);
    spinor_fraction operator/(const spinor &s1, const spinor_fraction &f2);
    spinor_fraction operator/(const spinor_fraction &f1, const spinor_product &p2);
    spinor_fraction operator/(const spinor_product &p1, const spinor_fraction &f2);
    spinor_fraction operator/(const spinor_fraction &f1, const spinor_sum &sum2);
    spinor_fraction operator/(const spinor_sum &sum1, const spinor_fraction &f2);
    spinor_fraction operator/(const spinor_fraction &f1, const spinor_fraction &f2);
    spinor_fraction operator/(const double &c, const spinor_fraction &f);
    spinor_fraction operator/(const spinor_fraction &f, const double &c);
    
    bool operator<(const spinor& a, const spinor& b);
    bool operator<(const spinor_product &p1, const spinor_product &p2);
    bool operator<(const factor_struct &fs1, const factor_struct &fs2);
    bool operator==(const spinor& a, const spinor& b);
    bool operator!=(const spinor& a, const spinor& b);
    bool operator==(const std::vector<spinor> &vec1, const std::vector<spinor> &vec2);
    bool operator==(const spinor_product& p1, const spinor_product& p2);
    bool operator!=(const spinor_product& p1, const spinor_product& p2);
    bool operator==(const spinor_sum& sum1, const spinor_sum& sum2);
    bool operator!=(const spinor_sum& sum1, const spinor_sum& sum2);

    bool is_similar_spinor(const spinor s1, const spinor s2);
    bool is_similar_spinor_vector(const std::vector<spinor> &vec1, const std::vector<spinor> &vec2);
    bool is_similar(const spinor_product p1, const spinor_product p2);
    bool is_spinor_in_numerator(const spinor_product &p, const spinor &s);
    bool is_number(std::string string);

    void cancellation(spinor_sum &sum);
    void order(spinor_product &p);
    void reduce_product(spinor_product &p);

    // Factorization functions
    std::vector<factor_struct> find_common_spinors(const spinor_sum &sum);
    int gcd(int a, int b);
    std::vector<int> find_unique_gcds(const spinor_sum &sum);
    std::set<std::vector<int>> make_gcd_combinations(const int &gcd_number);
    bool is_perfect_product(const spinor_sum &sum, const int &highest_order, const spinor_sum &factor,spinor_sum &leftover_sum);
    bool is_factor(const spinor_sum & sum, const spinor_sum &factor, spinor_sum &leftover_sum);
    factorized_spinor factor(factorized_spinor &sum, std::vector<factor_struct> factors);
    void factor(spinor_fraction &f);
    void flatten_fraction(spinor_fraction &f);
    factorized_spinor factor(const spinor &s, const spinor_sum &sum);
    void reduce_fraction(spinor_fraction &f);

    // Struct for long division
    struct Polynomial {
        std::map<int, spinor_sum> terms;
        int highest_order = 0;
        spinor var;
        
        Polynomial(const spinor_sum &sum, const spinor &var) {
            this->var = var;
            for (spinor_product prod : sum.terms) {
                int order = 0;
                for(spinor s : prod.numerator){
                    if (is_similar_spinor(s,var)){
                        order += 1;
                    }
                }

                terms[order].terms.push_back(prod);

                if (order > highest_order){
                    highest_order = order;
                }
            }
        }

        Polynomial operator-(const spinor_sum &sum){
            Polynomial result = *this;
            for (spinor_product prod : sum.terms) {
                int order = 0;
                for(spinor s : prod.numerator){
                    if (is_similar_spinor(s,var)){
                        order += 1;
                    }
                }

                result.terms[order] = result.terms[order] - prod;
            }
            return result;
        }

        bool is_null(){
            for(int i = 0; i < highest_order; i++){
                if(!terms[i].is_null){
                    return false;
                }
            }
            return true;
        }
    };


    std::string display_spinor(const spinor& s);

    std::ostream& operator<<(std::ostream& os, const spinor& s);
    std::ostream& operator<<(std::ostream& os, const spinor_product& p);
    std::ostream& operator<<(std::ostream& os, const spinor_sum& sum);
    std::ostream& operator<<(std::ostream& os,const spinor_fraction& f);
    std::ostream& operator<<(std::ostream& os,const factorized_spinor& f);

    std::vector<std::string> tokenize(const std::string &input);
    spinor parse_spinor(const std::string &token);

}