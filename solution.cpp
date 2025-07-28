#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <algorithm>

// --- Third-party library for JSON parsing ---
#include "json.hpp" 

// --- A lightweight BigInt library included directly in the file ---
// Source: A common, simplified BigInt implementation for competitive programming
// Note: This is a basic implementation for the purpose of this assignment.
#include <vector>
#include <string>
#include <ostream>
#include <istream>

struct BigInt {
    std::string a; 
    int sign;

    BigInt() : sign(1) {}

    BigInt(long long v) {
        *this = v;
    }

    BigInt(const std::string &s) {
        *this = s;
    }

    void operator=(long long v) {
        sign = (v < 0) ? -1 : 1;
        if (v < 0) v = -v;
        a.clear();
        if (v == 0) {
            a = "0";
            return;
        }
        while (v > 0) {
            a.push_back(v % 10 + '0');
            v /= 10;
        }
        std::reverse(a.begin(), a.end());
    }

    BigInt(const char* s) { *this = std::string(s); }

    void operator=(const std::string &s) {
        if (s.empty()) {
            sign = 1;
            a = "0";
            return;
        }
        sign = (s[0] == '-') ? -1 : 1;
        a = (sign == -1) ? s.substr(1) : s;
        if (a.empty() || a.find_first_not_of("0123456789") != std::string::npos) {
            throw std::invalid_argument("Invalid BigInt string");
        }
        if (a == "0") sign = 1;
    }

    friend BigInt operator+(BigInt a, BigInt b);
    friend BigInt operator-(BigInt a, BigInt b);
    friend BigInt operator*(BigInt a, BigInt b);
    friend BigInt operator/(BigInt a, BigInt b);

    friend std::ostream& operator<<(std::ostream &out, const BigInt &v) {
        if (v.sign == -1) out << '-';
        out << v.a;
        return out;
    }
    
    bool isZero() const { return a == "0"; }
};

// Basic implementations for +, -, * operators (simplified)
bool isSmaller(std::string str1, std::string str2) {
    int n1 = str1.length(), n2 = str2.length();
    if (n1 < n2) return true;
    if (n2 < n1) return false;
    for (int i = 0; i < n1; i++)
        if (str1[i] < str2[i]) return true;
        else if (str1[i] > str2[i]) return false;
    return false;
}

std::string findSum(std::string str1, std::string str2) {
    if (str1.length() > str2.length()) std::swap(str1, str2);
    std::string str = "";
    int n1 = str1.length(), n2 = str2.length();
    std::reverse(str1.begin(), str1.end());
    std::reverse(str2.begin(), str2.end());
    int carry = 0;
    for (int i=0; i<n1; i++) {
        int sum = ((str1[i]-'0')+(str2[i]-'0')+carry);
        str.push_back(sum%10 + '0');
        carry = sum/10;
    }
    for (int i=n1; i<n2; i++) {
        int sum = ((str2[i]-'0')+carry);
        str.push_back(sum%10 + '0');
        carry = sum/10;
    }
    if (carry) str.push_back(carry+'0');
    std::reverse(str.begin(), str.end());
    return str;
}

std::string findDiff(std::string str1, std::string str2) {
    std::string str = "";
    int n1 = str1.length(), n2 = str2.length();
    std::reverse(str1.begin(), str1.end());
    std::reverse(str2.begin(), str2.end());
    int carry = 0;
    for (int i=0; i<n2; i++) {
        int sub = ((str1[i]-'0')-(str2[i]-'0')-carry);
        if (sub < 0) {
            sub = sub + 10;
            carry = 1;
        } else carry = 0;
        str.push_back(sub + '0');
    }
    for (int i=n2; i<n1; i++) {
        int sub = ((str1[i]-'0') - carry);
        if (sub < 0) {
            sub = sub + 10;
            carry = 1;
        } else carry = 0;
        str.push_back(sub + '0');
    }
    std::reverse(str.begin(), str.end());
    // remove leading zeros
    size_t first_digit = str.find_first_not_of('0');
    if (std::string::npos != first_digit) return str.substr(first_digit);
    return "0";
}


BigInt operator+(BigInt a, BigInt b) {
    BigInt res;
    if (a.sign == b.sign) {
        res.sign = a.sign;
        res.a = findSum(a.a, b.a);
    } else {
        if (isSmaller(a.a, b.a)) std::swap(a, b);
        res.a = findDiff(a.a, b.a);
        res.sign = a.sign;
    }
    if (res.a == "0") res.sign = 1;
    return res;
}

BigInt operator-(BigInt a, BigInt b) {
    b.sign *= -1;
    return a + b;
}

BigInt operator*(BigInt a, BigInt b) {
    // FIX: Use 0LL to specify a long long literal, removing ambiguity.
    if (a.isZero() || b.isZero()) return BigInt(0LL);
    
    std::string s1 = a.a;
    std::string s2 = b.a;
    std::reverse(s1.begin(), s1.end());
    std::reverse(s2.begin(), s2.end());

    std::vector<int> m(s1.length() + s2.length(), 0);
    for (size_t i = 0; i < s1.length(); i++) {
        for (size_t j = 0; j < s2.length(); j++) {
            m[i + j] += (s1[i] - '0') * (s2[j] - '0');
        }
    }

    std::string product = "";
    for (size_t i = 0; i < m.size(); i++) {
        int digit = m[i] % 10;
        int carry = m[i] / 10;
        if (i + 1 < m.size()) {
            m[i + 1] += carry;
        }
        product = std::to_string(digit) + product;
    }
    
    size_t first_digit = product.find_first_not_of('0');
    if (std::string::npos != first_digit) {
       product = product.substr(first_digit);
    } else {
       product = "0";
    }

    BigInt res;
    res.a = product;
    res.sign = a.sign * b.sign;
    return res;
}


BigInt operator/(BigInt a, BigInt b) {
    if (b.isZero()) throw std::runtime_error("Division by zero");
    
    // FIX: Use 0LL to specify a long long literal.
    if (isSmaller(a.a, b.a)) return BigInt(0LL);
    
    // FIX: Explicitly cast the result to long long.
    if (a.a == b.a) return BigInt((long long)a.sign * b.sign);

    BigInt res;
    res.sign = a.sign * b.sign;
    
    std::string temp;
    std::string quotient;
    for (char digit : a.a) {
        temp += digit;
        int count = 0;
        BigInt temp_big(temp);
        BigInt b_abs = b; b_abs.sign = 1;
        while (!isSmaller(temp_big.a, b_abs.a)) {
            temp_big = temp_big - b_abs;
            count++;
        }
        quotient += std::to_string(count);
        temp = temp_big.a;
    }
    size_t first_digit = quotient.find_first_not_of('0');
    if (std::string::npos != first_digit) {
       res.a = quotient.substr(first_digit);
    } else {
       res.a = "0";
    }
    return res;
}


// Use nlohmann::json for convenience
using json = nlohmann::json;

// Function to convert value from a given base to a BigInt
BigInt string_to_bigint(const std::string& s, int base) {
    // FIX: Initialize with an explicit long long literal 0LL.
    BigInt res(0LL);
    BigInt big_base((long long)base);

    for (char c : s) {
        int digit;
        if (c >= '0' && c <= '9') {
            digit = c - '0';
        } else {
            digit = c - 'a' + 10;
        }
        res = res * big_base + BigInt((long long)digit);
    }
    return res;
}

// Function to solve for the secret 'c' given a JSON file
BigInt find_secret(const std::string& filename) {
    // 1. Read and parse the JSON file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    json data = json::parse(file);

    int k = data["keys"]["k"];
    std::vector<BigInt> x_coords;
    std::vector<BigInt> y_coords;

    // 2. Decode the points (x, y)
    int count = 0;
    for (auto& [key, val] : data.items()) {
        if (key == "keys") continue;

        if (count < k) {
            // x is the key
            x_coords.push_back(BigInt(stoll(key)));

            // y is the value, decoded from its base
            int base = std::stoi(val["base"].get<std::string>());
            std::string value_str = val["value"].get<std::string>();
            y_coords.push_back(string_to_bigint(value_str, base));
            count++;
        }
    }

    // 3. Find the Secret (c) using Lagrange Interpolation at x=0
    // The constant term c is f(0).
    // f(0) = sum_{j=0}^{k-1} y_j * L_j(0)
    // L_j(0) = product_{i=0, i!=j}^{k-1} x_i / (x_i - x_j)

    // FIX: Initialize with an explicit long long literal 0LL.
    BigInt secret(0LL);

    for (int j = 0; j < k; ++j) {
        // FIX: Initialize with an explicit long long literal 1LL.
        BigInt numerator(1LL);
        BigInt denominator(1LL);

        for (int i = 0; i < k; ++i) {
            if (i == j) continue;
            numerator = numerator * x_coords[i];
            denominator = denominator * (x_coords[i] - x_coords[j]);
        }
        
        // The Lagrange basis polynomial evaluated at 0
        BigInt lambda_j = numerator / denominator;
        
        // Add this part of the sum to the total secret
        secret = secret + (y_coords[j] * lambda_j);
    }

    return secret;
}


int main() {
    try {
        std::string filename;
        std::cout << "Enter JSON filename (e.g., testcase1.json): ";
        std::cin >> filename;

        BigInt secret = find_secret(filename);
        std::cout << "Secret: " << secret << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
