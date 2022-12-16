#pragma once
#include <algorithm>
#include <compare>
#include <complex>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <vector>

class BigInteger {
 private:
  std::vector<int> digits_;
  bool is_positive_;

  static void fft(std::vector<std::complex<double>>&, bool);
  static size_t reverseBits(size_t, size_t);
  static size_t closestPower(size_t);
  static size_t getBinPower(size_t);
  std::vector<std::complex<double>> convertDigits() const;
  void fixDigits();
  BigInteger stupidMult(int) const;

 public:
  static const int BASE = 1'00;
  static const int DEC_POW = 2;

  BigInteger();
  BigInteger(int number);
  BigInteger(const BigInteger&) = default;
  explicit BigInteger(const std::string&);

  BigInteger operator+() const;
  BigInteger operator-() const;
  BigInteger& operator++();
  BigInteger operator++(int);
  BigInteger& operator--();
  BigInteger operator--(int);
  friend void swap(BigInteger&, BigInteger&);
  BigInteger& operator=(const BigInteger&) = default;
  BigInteger& operator+=(const BigInteger&);
  BigInteger& operator-=(const BigInteger&);
  BigInteger& operator*=(const BigInteger&);
  //BigInteger& operator/=(const BigInteger&);
  friend BigInteger& operator/=(BigInteger&, const BigInteger&);
  BigInteger& operator%=(const BigInteger&);

  std::string toString(bool) const;
  const std::vector<int>& digits() const;
  size_t digitsSize() const;
  void changeSign();
  bool isPositive() const;
  bool isZero() const;
  std::strong_ordering absCmp(const BigInteger&) const;
  static BigInteger gcd(const BigInteger&, const BigInteger&);
  void shiftBase(int);

  explicit operator bool() const;
};

class Rational {
 private:
  BigInteger numerator_, denominator_;

  void fixDividers();
  bool isZero() const;

 public:
  Rational();
  Rational(int, int);
  Rational(const BigInteger&, const BigInteger&);
  Rational(const Rational&) = default;

  void changeSign();
  bool isPositive() const;
  const BigInteger& numerator() const;
  const BigInteger& denominator() const;

  Rational operator+();
  Rational operator-();
  friend void swap(Rational&, Rational&);
  Rational& operator=(const Rational&) = default;
  Rational& operator+=(const Rational&);
  Rational& operator-=(const Rational&);
  Rational& operator*=(const Rational&);
  Rational& operator/=(const Rational&);

  std::string toString() const;
  std::string asDecimal(size_t) const;
  explicit operator double() const;
};

BigInteger operator""_bi(unsigned long long number) {
  return BigInteger(std::to_string(number));
}

BigInteger operator""_bi(const char* c_str) {
  std::string s = c_str;
  return BigInteger(s);
}

BigInteger operator""_bi(const char* c_str, size_t) {
  return BigInteger(c_str);
}

size_t BigInteger::getBinPower(size_t number) {
  for (size_t i = 0; i < INT32_WIDTH; ++i) {
    if (static_cast<size_t>(1 << i) >= number) return i;
  }
  return INT32_WIDTH;
}

size_t BigInteger::reverseBits(size_t number, size_t size) {
  for (size_t i = 0; i < size / 2; ++i) {
    size_t get = static_cast<size_t>((1 << i) + (1 << (size - i - 1)));
    get &= number;
    get = ((get << (size - 2 * i - 1)) & (1 << (size - i - 1))) +
          ((get >> (size - 2 * i - 1)) & (1 << i));
    number ^= number & static_cast<size_t>(((1 << i) + (1 << (size - i - 1))));
    number |= get;
  }
  return number;
}

size_t BigInteger::closestPower(size_t number) {
  for (size_t i = 0; i < INT32_WIDTH; ++i) {
    if (static_cast<size_t>(1 << i) >= number) return (1 << i);
  }
  return INT32_WIDTH;
}

std::vector<std::complex<double>> BigInteger::convertDigits() const {
  std::vector<std::complex<double>> complex(digitsSize());
  for (size_t index = 0; index < digitsSize(); ++index) {
    complex[index] = digits_[index];
  }
  size_t closest_power = BigInteger::closestPower(complex.size());
  while (complex.size() != closest_power) {
    complex.push_back(0.0);
  }
  return complex;
}

void BigInteger::fft(std::vector<std::complex<double>>& coeffs,
                     bool reverse = false) {
  size_t size = coeffs.size();

  std::complex<double> quark =
      std::polar(1.0, 2 * M_PI / static_cast<double>(size));
  if (reverse) {
    quark = std::complex<double>(1.0) / quark;
  }
  if (size == 1) {
    return;
  }

  for (size_t i = 0; i < size; ++i) {
    size_t reversed_i = reverseBits(i, getBinPower(size));
    if (i < reversed_i) std::swap(coeffs[i], coeffs[reversed_i]);
  }
  size_t end = size;
  for (size_t step = 2; step <= size; step *= 2) {
    std::complex<double> current_quark = quark;
    for (size_t i = size; i > step; i /= 2) current_quark *= current_quark;
    for (size_t start = 0; start != end; start += step) {
      size_t middle = start + step / 2;
      std::complex<double> quark_degree = 1.0;
      for (size_t left = start, right = middle; left != middle;
           left++, right++) {
        std::complex<double> first = coeffs[left];
        std::complex<double> second = coeffs[right] * quark_degree;
        coeffs[left] = first + second;
        coeffs[right] = first - second;
        quark_degree *= current_quark;
      }
    }
  }
  if (reverse) {
    for (size_t i = 0; i < size; ++i) coeffs[i] /= static_cast<double>(size);
  }
}

BigInteger::BigInteger()
    : digits_(std::vector<int>(1, 0)), is_positive_(true) {}
BigInteger::BigInteger(int number)
    : digits_(std::vector<int>(1, std::abs(number))),
      is_positive_(number >= 0) {
  fixDigits();
}

BigInteger::BigInteger(const std::string& inp_str)
    : digits_(std::vector<int>()), is_positive_(true) {
  std::string buf(DEC_POW, '0');
  std::string str = inp_str;
  buf[DEC_POW] = 0;
  size_t index = 0;

  is_positive_ = (str[index] != '-');
  index = (str[index] == '-') || (str[index] == '+');
  std::reverse(str.begin() + static_cast<int>(index), str.end());

  while (index < str.size()) {
    size_t buf_index = 0;
    std::fill(buf.begin(), buf.end(), '0');
    while (index < str.size() && buf_index < buf.size()) {
      buf[buf_index] = str[index];
      ++index;
      ++buf_index;
    }
    buf[buf_index] = 0;
    std::reverse(buf.begin(), buf.begin() + static_cast<int>(buf_index));
    digits_.push_back(std::stoi(buf));
  }
  fixDigits();
}

const std::vector<int>& BigInteger::digits() const { return digits_; }

size_t BigInteger::digitsSize() const { return digits_.size(); }

void BigInteger::changeSign() {
  is_positive_ = !is_positive_;
  if (isZero()) is_positive_ = true;
}

bool BigInteger::isPositive() const { return is_positive_; }

bool BigInteger::isZero() const {
  return digits_.size() == 1 && digits_[0] == 0;
}

BigInteger abs(const BigInteger& arg) {
  BigInteger tmp = arg;
  if (!arg.isPositive()) tmp.changeSign();
  return tmp;
}

void BigInteger::shiftBase(int shift_size) {
  if (shift_size >= 0) {
    BigInteger shift;
    shift.digits_.resize(static_cast<size_t>(shift_size + 1), 0);
    *(shift.digits_.end() - 1) = 1;
    *this *= shift;
    return;
  }
  if (std::abs(shift_size) >= static_cast<int>(digitsSize())) {
    digits_.resize(1, 1);
    return;
  }
  std::reverse(digits_.begin(), digits_.end());
  digits_.resize(digitsSize() + static_cast<size_t>(shift_size));
  std::reverse(digits_.begin(), digits_.end());
}

BigInteger::operator bool() const { return !isZero(); }

void BigInteger::fixDigits() {
  while (digitsSize() > 1 && digits_.back() == 0) {
    digits_.pop_back();
  }

  if (*(digits_.rbegin()) < 0) {
    is_positive_ = !is_positive_;
    for (size_t index = 0; index < digitsSize(); ++index) {
      digits_[index] = -digits_[index];
    }
  }

  for (size_t index = digitsSize() - 1; index > 0; --index) {
    if (digits_[index - 1] <= 0) {
      int add = std::abs(digits_[index - 1]) / BASE + 1;
      digits_[index] -= add;
      digits_[index - 1] += add * BASE;
    }
  }

  int go = 0, stay = 0;
  for (size_t index = 0; index < digitsSize(); ++index) {
    digits_[index] += go;
    stay = digits_[index] % BASE;
    go = digits_[index] / BASE;
    digits_[index] = stay;
  }
  while (go != 0) {
    digits_.push_back(go % BASE);
    go /= BASE;
  }

  while (digitsSize() > 1 && digits_.back() == 0) {
    digits_.pop_back();
  }

  if (digits_.size() == 1 && digits_[0] == 0) is_positive_ = true;
}

std::strong_ordering BigInteger::absCmp(const BigInteger& other) const {
  auto size_ordering = digitsSize() <=> other.digitsSize();
  if (size_ordering != std::strong_ordering::equal) return size_ordering;
  for (size_t index = digitsSize(); index > 0; --index) {
    if (digits()[index - 1] > other.digits()[index - 1])
      return std::strong_ordering::greater;
    if (digits()[index - 1] < other.digits()[index - 1])
      return std::strong_ordering::less;
  }
  return std::strong_ordering::equal;
}

void swap(BigInteger& lhs, BigInteger& rhs) {
  std::swap(lhs.is_positive_, rhs.is_positive_);
  std::swap(lhs.digits_, rhs.digits_);
}

bool operator==(const BigInteger& lhs, const BigInteger& rhs) {
  if (&lhs == &rhs) return true;
  if (lhs.isZero() && rhs.isZero()) return true;
  if (lhs.digitsSize() != rhs.digitsSize()) return false;
  for (size_t index = 0; index < lhs.digitsSize(); ++index) {
    if (lhs.digits()[index] != rhs.digits()[index]) return false;
  }
  if (lhs.isPositive() != rhs.isPositive()) return false;
  return true;
}

auto operator<=>(const BigInteger& lhs, const BigInteger& rhs) {
  if (lhs.isPositive() ^ rhs.isPositive())
    return lhs.isPositive() <=> rhs.isPositive();
  if (lhs == rhs) return std::strong_ordering::equal;
  return (lhs.isPositive() && rhs.isPositive()) ? lhs.absCmp(rhs)
                                                : 0 <=> lhs.absCmp(rhs);
}

BigInteger& BigInteger::operator+=(const BigInteger& other) {
  for (size_t i = 0; i < other.digitsSize(); ++i) {
    if (digitsSize() == i) {
      digits_.push_back(0);
    }
    digits_[i] += (is_positive_ == other.is_positive_ ? other.digits_[i]
                                                      : -other.digits_[i]);
  }
  fixDigits();
  return *this;
}

BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  return tmp += rhs;
}

BigInteger& BigInteger::operator-=(const BigInteger& other) {
  changeSign();
  *this += other;
  changeSign();
  return *this;
}

BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  return tmp -= rhs;
}

BigInteger BigInteger::operator+() const { return *this; }
BigInteger BigInteger::operator-() const {
  BigInteger tmp = *this;
  tmp.changeSign();
  return tmp;
}

BigInteger& BigInteger::operator++() {
  digits_[0] += (is_positive_ ? 1 : -1);
  if (digits_[0] >= static_cast<int>(BASE)) fixDigits();
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger tmp = *this;
  ++(*this);
  return tmp;
}

BigInteger& BigInteger::operator--() {
  digits_[0] -= (is_positive_ ? 1 : -1);
  if (digits_[0] < 0) fixDigits();
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger tmp = *this;
  --(*this);
  return tmp;
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
  std::vector<std::complex<double>> left = convertDigits(),
                                    right = other.convertDigits();
  left.resize(std::max(left.size(), right.size()) * 2, 0.0);
  right.resize(left.size(), 0.0);
  BigInteger::fft(left);
  BigInteger::fft(right);
  for (size_t index = 0; index < left.size(); ++index) {
    left[index] *= right[index];
  }
  BigInteger::fft(left, true);
  digits_.resize(left.size());
  long long digit = 0;
  const double MAGIC_NUMBER = 0.5;
  for (size_t index = 0; index < left.size(); ++index) {
    digit +=
        static_cast<long long>(std::floor(left[index].real() + MAGIC_NUMBER));
    digits_[index] = static_cast<int>(digit % BigInteger::BASE);
    digit /= BigInteger::BASE;
  }
  if (digit) {
    digits_.push_back(static_cast<int>(digit));
  }

  is_positive_ = !(is_positive_ ^ other.is_positive_);
  fixDigits();
  return *this;
}

BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  tmp *= rhs;
  return tmp;
}

std::string BigInteger::toString(bool leading_zeros = false) const {
  std::stringstream ss;
  if (!is_positive_) ss << "-";
  if (!leading_zeros) {
    ss << digits_[digits_.size() - 1];
  } else {
    ss << std::setw(DEC_POW) << std::setfill('0')
       << digits_[digits_.size() - 1];
  }
  if (digits_.size() > 1) {
    for (size_t index = digits_.size() - 1; index > 0; --index) {
      ss << std::setw(DEC_POW) << std::setfill('0') << digits_[index - 1];
    }
  }
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, const BigInteger& rhs) {
  os << rhs.toString();
  return os;
}

std::istream& operator>>(std::istream& is, BigInteger& rhs) {
  std::string str;
  is >> str;
  rhs = BigInteger(str);
  return is;
}

BigInteger BigInteger::stupidMult(int number) const {
  BigInteger tmp = *this;
  for (size_t i = 0; i < tmp.digitsSize(); ++i) tmp.digits_[i] *= number;
  tmp.fixDigits();
  return tmp;
}

BigInteger& operator/=(BigInteger& lhs, const BigInteger& rhs) {
  // // std::cerr << "/=";
  if (&lhs == &rhs) return lhs = 1;
  if (lhs.absCmp(rhs) == std::strong_ordering::less) return lhs = 0;
  if (rhs.isZero()) {
    throw "BigInteger: Division by zero.";
  }
  BigInteger calc;
  BigInteger ans;
  BigInteger shift;
  bool l_sign = lhs.is_positive_, r_sign = rhs.is_positive_;
  ans.digits_.resize(0);
  shift.digits_.resize(lhs.digitsSize() - rhs.digitsSize() + 1, 0);

  while (lhs.absCmp(rhs) >= 0) {
    calc.digits_.resize(rhs.digitsSize());
    std::copy(lhs.digits_.begin() + static_cast<int>(lhs.digitsSize()) -
                  static_cast<int>(calc.digitsSize()),
              lhs.digits_.end(), calc.digits_.begin());
    while (calc.absCmp(rhs) == std::strong_ordering::less) {
      calc.digits_.insert(
          calc.digits_.begin(),
          lhs.digits_[lhs.digitsSize() - calc.digitsSize() - 1]);
    }

    for (size_t i = 0;
         i + lhs.digitsSize() + 2 < shift.digitsSize() + calc.digitsSize(); ++i)
      ans.digits_.push_back(0);

    int left = 0, right = BigInteger::BASE;
    while (left + 1 < right) {
      int mid = (left + right) / 2;
      if (calc.absCmp(abs(rhs.stupidMult(mid))) == std::strong_ordering::less) {
        right = mid;
      } else {
        left = mid;
      }
    }
    int digit = left;

    shift.digits_.resize(lhs.digitsSize() - calc.digitsSize() + 1, 0);
    BigInteger rhs_c;
    rhs_c.digits_.resize(shift.digitsSize() - 1 + rhs.digitsSize(), 0);
    std::copy(rhs.digits_.begin(), rhs.digits_.end(),
              rhs_c.digits_.begin() + static_cast<int>(shift.digitsSize() - 1));
    *(shift.digits_.end() - 1) = 1;

    lhs -=
        (l_sign ? abs(rhs_c.stupidMult(digit)) : -abs(rhs_c.stupidMult(digit)));
    ans.digits_.push_back(digit);
  }
  std::reverse(ans.digits_.begin(), ans.digits_.end());
  ans *= shift;
  lhs.digits_ = ans.digits_;
  lhs.is_positive_ = !(l_sign ^ r_sign);
  lhs.fixDigits();
  return lhs;
}

BigInteger operator/(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  tmp /= rhs;
  return tmp;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
  return *this = *this - *this / other * other;
}

BigInteger operator%(const BigInteger& lhs, const BigInteger& rhs) {
  BigInteger tmp = lhs;
  return tmp %= rhs;
}

Rational::Rational() : numerator_(0ll), denominator_(1ll) {}

Rational::Rational(int numerator, int denominator = 1ll)
    : numerator_(numerator), denominator_(denominator) {
  fixDividers();
}

Rational::Rational(const BigInteger& numerator,
                   const BigInteger& denominator = BigInteger(1ll))
    : numerator_(numerator), denominator_(denominator) {
  fixDividers();
}

void Rational::changeSign() { numerator_.changeSign(); }

bool Rational::isPositive() const { return numerator_.isPositive(); }

const BigInteger& Rational::numerator() const { return numerator_; }

const BigInteger& Rational::denominator() const { return denominator_; }

BigInteger BigInteger::gcd(const BigInteger& lhs, const BigInteger& rhs) {
  if (rhs.isZero()) return abs(lhs);
  return gcd(abs(rhs), abs(lhs) % abs(rhs));
}

void Rational::fixDividers() {
  if (denominator_.isZero()) throw "Rational: Division by zero.";
  if (!denominator_.isPositive()) {
    numerator_.changeSign();
    denominator_.changeSign();
  }
  BigInteger gcd = BigInteger::gcd(numerator_, denominator_);
  numerator_ /= gcd;
  denominator_ /= gcd;
}

Rational Rational::operator+() { return *this; }
Rational Rational::operator-() {
  Rational tmp = *this;
  tmp.changeSign();
  return tmp;
}

void swap(Rational& lhs, Rational& rhs) {
  swap(lhs.numerator_, rhs.numerator_);
  swap(lhs.denominator_, rhs.denominator_);
}

Rational& Rational::operator+=(const Rational& other) {
  numerator_ =
      numerator_ * other.denominator_ + other.numerator_ * denominator_;
  denominator_ *= other.denominator_;
  fixDividers();
  return *this;
}

Rational operator+(const Rational& lhs, const Rational& rhs) {
  Rational tmp = lhs;
  tmp += rhs;
  return tmp;
}

Rational& Rational::operator-=(const Rational& other) {
  numerator_ =
      numerator_ * other.denominator_ - other.numerator_ * denominator_;
  denominator_ *= other.denominator_;
  fixDividers();
  return *this;
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
  Rational tmp = lhs;
  return tmp -= rhs;
}

Rational& Rational::operator*=(const Rational& other) {
  numerator_ *= other.numerator_;
  denominator_ *= other.denominator_;
  fixDividers();
  return *this;
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
  Rational tmp = lhs;
  return tmp *= rhs;
}

Rational& Rational::operator/=(const Rational& other) {
  if (this == &other) return *this = 1;
  numerator_ *= other.denominator_;
  denominator_ *= other.numerator_;
  fixDividers();
  return *this;
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
  Rational tmp = lhs;
  return tmp /= rhs;
}

bool operator==(const Rational& lhs, const Rational& rhs) {
  return lhs.numerator() == rhs.numerator() &&
         lhs.denominator() == rhs.denominator();
}

auto operator<=>(const Rational& lhs, const Rational& rhs) {
  return (lhs.numerator() * rhs.denominator()) <=>
         (rhs.numerator() * lhs.denominator());
}

std::string Rational::toString() const {
  return numerator_.toString() +
         (denominator_ != BigInteger(1ll) ? "/" + denominator_.toString() : "");
}

std::string Rational::asDecimal(size_t precision) const {
  if (numerator_.isZero()) return "0." + std::string(precision, '0');
  std::string ret = "";
  BigInteger part = abs(numerator_) / denominator_;
  if (!numerator_.isPositive()) ret += "-";
  ret += part.toString();

  part = abs(numerator_) % denominator_;
  part.shiftBase(2 * static_cast<int>(precision) / BigInteger::DEC_POW);
  part /= denominator_;
  std::string parts = part.toString(true);
  std::string fract(
      2 * precision / BigInteger::DEC_POW * BigInteger::DEC_POW - parts.size(),
      '0');
  fract += parts;
  fract.resize(precision + 1, '0');
  if (fract[fract.size() - 1] >= '5') fract[fract.size() - 2] += 1;
  fract.pop_back();
  ret += "." + fract;
  return ret;
}

Rational::operator double() const {
  const size_t precision = 40;
  return std::stod(asDecimal(precision));
}
