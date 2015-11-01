/*******************************************************************************
  SIDIM - SI-dimension

  Copyright (c) 2015, Manuel Freiberger
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  - Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
  - Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

#ifndef SIDIM_HPP
#define SIDIM_HPP

//#define SIDIM_USE_WEOS

#ifdef SIDIM_USE_WEOS

#include <weos/ratio.hpp>
#include <weos/type_traits.hpp>
#define SIDIM_STD weos

#else

#include <ratio>
#include <type_traits>
#define SIDIM_STD std

#endif // SIDIM_USE_WEOS

#include <cstdint>
#include <limits>


namespace sidim
{

//! A distance.
//!
//! The \p distance type is used to quantify differences between two positions
//! in 1D-space. The accuracy is determined via the template parameter \p TRep
//! and the units via \p TInterval. \p TRep must be a signed and can be an
//! integer or a floating point type. The interval type is a compile-time ratio
//! and specifies the interval between two distance values in fractions of
//! a meter. For example,
//! \code
//! distance<int32_t, std::milli>
//! \endcode
//! is a distance type representing millimeters whose values are represented
//! by 32-bit integers.
//!
//! Distances are implicitly convertible to other types iff the conversion
//! is lossless. For example, a <tt>distance<int, std::micro></tt> can be
//! constructed from a <tt>distance<int, std::milli></tt> because the former
//! is a more accurate representation. To convert in the other direction,
//! make use of \p distance_cast.
template <typename TRep, typename TInterval = SIDIM_STD::ratio<1>>
class distance;

//! A position in 1D-space.
//!
//! The position is a position in 1D-space. The \p TTag template parameter
//! allows to create different position types which are unrelated to each
//! other, e.g. to quantify positions independently along the x- and y-axis.
//! The \p TDistance template parameter determines the units which are used
//! for the representation of a position.
//! The following lines define types to measure a position in x-direction
//! in micrometers and in y-direction in meters.
//! \code
//! struct x_tag;
//! using x_position_t = position<x_tag, distance<int, std::micro>>;
//! struct y_tag;
//! using y_position_t = position<y_tag, distance<int, std:ratio<1>>>;
//! \endcode
//!
//! All positions are measured from an origin.
template <typename TTag, typename TDistance>
class position;

namespace sidim_detail
{

template <typename TFromDistance, typename TToDistance>
struct distance_caster;

// -----------------------------------------------------------------------------
// is_distance
// -----------------------------------------------------------------------------

// A trait to determine if a type is a distance.
template <typename TType>
struct is_distance : public SIDIM_STD::false_type {};

template <typename TRep, typename TInterval>
struct is_distance<distance<TRep, TInterval>> : public SIDIM_STD::true_type {};

template <typename TRep, typename TInterval>
struct is_distance<const distance<TRep, TInterval>> : public SIDIM_STD::true_type {};

template <typename TRep, typename TInterval>
struct is_distance<volatile distance<TRep, TInterval>> : public SIDIM_STD::true_type {};

template <typename TRep, typename TInterval>
struct is_distance<const volatile distance<TRep, TInterval>> : public SIDIM_STD::true_type {};

// -----------------------------------------------------------------------------
// Compile-time Greatest Common Divisor & Least Common Multiple
// -----------------------------------------------------------------------------

// EUCLID(a, b)
//   if b = 0
//     then return a
//   else return EUCLID(b, a mod b)
template <intmax_t A, intmax_t B>
struct static_gcd_impl
{
    static const intmax_t value = static_gcd_impl<B, A % B>::value;
};

template <intmax_t A>
struct static_gcd_impl<A, 0>
{
    static const intmax_t value = A;
};

template <>
struct static_gcd_impl<0, 0>
{
    static const intmax_t value = 0;
};

template <intmax_t A, intmax_t B>
struct static_gcd : public SIDIM_STD::integral_constant<intmax_t, static_gcd_impl<A, B>::value>
{
};


template <intmax_t A, intmax_t B>
struct static_lcm_impl
{
    static const intmax_t value = A / static_gcd_impl<A, B>::value * B;
};

template <>
struct static_lcm_impl<0, 0>
{
    static const intmax_t value = 0;
};

template <intmax_t A, intmax_t B>
struct static_lcm : public SIDIM_STD::integral_constant<intmax_t, static_lcm_impl<A, B>::value>
{
};

// -----------------------------------------------------------------------------
// ratio_gcd
// -----------------------------------------------------------------------------

template <typename R1, typename R2>
struct ratio_gcd : SIDIM_STD::ratio<static_gcd<R1::num, R2::num>::value,
                                   static_lcm<R1::den, R2::den>::value>::type
{
};

// -----------------------------------------------------------------------------
// checked_division
// -----------------------------------------------------------------------------

template <typename R1, typename R2>
class checked_division
{
    // Divide the numerators by there GCD.
    static const intmax_t gcd_num = static_gcd<R1::num, R2::num>::value;
    static const intmax_t num1 = R1::num / gcd_num;
    static const intmax_t num2 = R2::num / gcd_num;
    // Do the same with the denominators.
    static const intmax_t gcd_den = static_gcd<R1::den, R2::den>::value;
    static const intmax_t den1 = R1::den / gcd_den;
    static const intmax_t den2 = R2::den / gcd_den;

    // Assume sizeof(intmax_t) == 2 and CHAR_BIT == 8:
    // We want max = 0x7fff.
    // Set the highest bit and add 1:
    //   1 << (16 - 1) + 1 = 0x4000 + 1 = 0x4001
    // Negation by two's complement is bit-wise negation plus 1:
    //   -0x4001 = ~0x4001 + 1 = 0x7ffe + 1 = 0x7fff
    // This is what we want.
    static const intmax_t max = -((intmax_t(1) << (sizeof(intmax_t) * 8 - 1)) + 1);

public:
    static const bool overflow = (num1 > max / den2) || (num2 > max / den1);

    // (n1 / d1) / (n2 / d2) = (n1 * d2) / (d1 * n2)
    typedef SIDIM_STD::ratio<num1 * den2, den1 * num2> type;
};

} // namespace sidim_detail
} // namespace sidim

// ----=====================================================================----
//     Specialization of common_type for sidim::distance<>
// ----=====================================================================----

namespace SIDIM_STD
{

template <typename TRep1, typename TInterval1,
          typename TRep2, typename TInterval2>
struct common_type<sidim::distance<TRep1, TInterval1>,
                   sidim::distance<TRep2, TInterval2>>
{
    using type = sidim::distance<
                     typename common_type<TRep1, TRep2>::type,
                     typename sidim::sidim_detail::ratio_gcd<TInterval1, TInterval2>::type>;
};

} // namespace SIDIM_STD

namespace sidim
{

// ----=====================================================================----
//     treat_as_floating_point
// ----=====================================================================----

template <typename TRep>
struct treat_as_floating_point : public SIDIM_STD::is_floating_point<TRep> {};

// ----=====================================================================----
//     distance_values
// ----=====================================================================----

template <typename TRep>
struct distance_values
{
    static constexpr
    TRep zero()
    {
        return TRep(0);
    }

    static constexpr
    TRep min()
    {
        return std::numeric_limits<TRep>::lowest();
    }

    static constexpr
    TRep max()
    {
        return std::numeric_limits<TRep>::max();
    }
};

// ----=====================================================================----
//     distance
// ----=====================================================================----

template <typename TToDistance, typename TRep, typename TInterval>
inline constexpr
typename SIDIM_STD::enable_if<sidim_detail::is_distance<TToDistance>::value,
                             TToDistance>::type
distance_cast(const distance<TRep, TInterval>& d);

template <typename TRep, typename TInterval>
class distance
{
    static_assert(!sidim_detail::is_distance<TRep>::value,
                  "The representation cannot be a distance");
    static_assert(TInterval::num > 0, "The interval must be positive");

public:
    //! The type used for representing the number of intervals.
    typedef TRep rep;
    //! The interval between two values specified as a fraction of a meter.
    typedef TInterval interval;

    //! Creates a distance.
    constexpr
    distance() = default;

    distance(const distance& other) = default;

    template <typename TRep2,
              typename = typename SIDIM_STD::enable_if<
                  SIDIM_STD::is_convertible<TRep2, rep>::value
                  && (treat_as_floating_point<rep>::value ||
                      !treat_as_floating_point<TRep2>::value)
              >::type>
    constexpr explicit
    distance(const TRep2& count)
        : m_count(count)
    {
    }

    template <typename TRep2, typename TInterval2,
              typename = typename SIDIM_STD::enable_if<
                  !sidim_detail::checked_division<TInterval2, interval>::overflow
                  && (treat_as_floating_point<rep>::value
                      || (sidim_detail::checked_division<TInterval2, interval>::type::den == 1
                      && !treat_as_floating_point<TRep2>::value))
              >::type>
    constexpr
    distance(const distance<TRep2, TInterval2>& other)
        : m_count(distance_cast<distance>(other).count())
    {
    }

    distance& operator=(const distance& other) = default;

    constexpr
    rep count() const
    {
        return m_count;
    }

    // Arithmetic operators.

    constexpr
    distance operator+() const
    {
        return *this;
    }

    constexpr
    distance operator-() const
    {
        return distance(-m_count);
    }

    distance& operator++()
    {
        ++m_count;
        return *this;
    }

    // \todo Other operators are missing

    distance operator++(int)
    {
        return distance(m_count++);
    }

    distance& operator--()
    {
        --m_count;
        return *this;
    }

    distance operator--(int)
    {
        return distance(m_count--);
    }

    distance& operator+=(const distance& other)
    {
        m_count += other.m_count;
        return *this;
    }

    distance& operator-=(const distance& other)
    {
        m_count -= other.m_count;
        return *this;
    }

    distance& operator*=(const rep& a)
    {
        m_count *= a;
        return *this;
    }

    distance& operator/=(const rep& a)
    {
        m_count /= a;
        return *this;
    }

    distance& operator%=(const rep& a)
    {
        m_count %= a;
        return *this;
    }

    distance& operator%=(const distance& other)
    {
        m_count %= other.count();
        return *this;
    }

    // Special values.

    static constexpr
    distance zero()
    {
        return distance(distance_values<rep>::zero());
    }

    static constexpr
    distance min()
    {
        return distance(distance_values<rep>::min());
    }

    static constexpr
    distance max()
    {
        return distance(distance_values<rep>::max());
    }

private:
    //! The number of intervals of which this distance consists.
    rep m_count;
};

// ----=====================================================================----
//     SI-constants
// ----=====================================================================----

typedef distance<std::int64_t, SIDIM_STD::nano>  nanometer;
typedef distance<std::int64_t, SIDIM_STD::micro> micrometer;
typedef distance<std::int64_t, SIDIM_STD::milli> millimeter;
typedef distance<std::int64_t>                  meter;

// ----=====================================================================----
//     distance comparisons
// ----=====================================================================----

namespace sidim_detail
{

// Compares two distances. General case.
template <typename TDistance1, typename TDistance2>
struct distance_equal
{
    static inline constexpr
    bool cmp(const TDistance1& d1, const TDistance2& d2)
    {
        using common_type = typename SIDIM_STD::common_type<TDistance1, TDistance2>::type;
        return common_type(d1).count() == common_type(d2).count();
    }
};

// Compares two distances. No cast needed when both distances are of the
// same type.
template <typename TDistance>
struct distance_equal<TDistance, TDistance>
{
    static inline constexpr
    bool cmp(const TDistance& d1, const TDistance& d2)
    {
        return d1.count() == d2.count();
    }
};

// Compares two distances. General case.
template <typename TDistance1, typename TDistance2>
struct distance_less
{
    static inline constexpr
    bool cmp(const TDistance1& d1, const TDistance2& d2)
    {
        using common_type = typename SIDIM_STD::common_type<TDistance1, TDistance2>::type;
        return common_type(d1).count() < common_type(d2).count();
    }
};

// Compares two distances. No cast needed when both distances are of the
// same type.
template <typename TDistance>
struct distance_less<TDistance, TDistance>
{
    static inline constexpr
    bool cmp(const TDistance& d1, const TDistance& d2)
    {
        return d1.count() < d2.count();
    }
};

} // namespace sidim_detail

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
bool operator==(const distance<TRep1, TInterval1>& d1, const distance<TRep2, TInterval2>& d2)
{
    return sidim_detail::distance_equal<distance<TRep1, TInterval1>,
                                       distance<TRep2, TInterval2> >::cmp(d1, d2);
}

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
bool operator!=(const distance<TRep1, TInterval1>& d1, const distance<TRep2, TInterval2>& d2)
{
    return !sidim_detail::distance_equal<distance<TRep1, TInterval1>,
                                        distance<TRep2, TInterval2> >::cmp(d1, d2);
}

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
bool operator<(const distance<TRep1, TInterval1>& d1, const distance<TRep2, TInterval2>& d2)
{
    return sidim_detail::distance_less<distance<TRep1, TInterval1>,
                                      distance<TRep2, TInterval2> >::cmp(d1, d2);
}

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
bool operator>(const distance<TRep1, TInterval1>& d1, const distance<TRep2, TInterval2>& d2)
{
    return sidim_detail::distance_less<distance<TRep2, TInterval2>,
                                      distance<TRep1, TInterval1> >::cmp(d2, d1);
}

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
bool operator<=(const distance<TRep1, TInterval1>& d1, const distance<TRep2, TInterval2>& d2)
{
    return !sidim_detail::distance_less<distance<TRep2, TInterval2>,
                                       distance<TRep1, TInterval1> >::cmp(d2, d1);
}

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
bool operator>=(const distance<TRep1, TInterval1>& d1, const distance<TRep2, TInterval2>& d2)
{
    return !sidim_detail::distance_less<distance<TRep1, TInterval1>,
                                       distance<TRep2, TInterval2> >::cmp(d1, d2);
}

// ----=====================================================================----
//     distance arithmetic
// ----=====================================================================----

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
typename SIDIM_STD::common_type<distance<TRep1, TInterval1>,
                               distance<TRep2, TInterval2> >::type
operator+(const distance<TRep1, TInterval1>& x,
          const distance<TRep2, TInterval2>& y)
{
    using common_type = typename SIDIM_STD::common_type<
            distance<TRep1, TInterval1>,
            distance<TRep2, TInterval2> >::type;
    return common_type(common_type(x).count() + common_type(y).count());
}

template <typename TRep1, typename TInterval1, typename TRep2, typename TInterval2>
inline constexpr
typename SIDIM_STD::common_type<distance<TRep1, TInterval1>,
                               distance<TRep2, TInterval2> >::type
operator-(const distance<TRep1, TInterval1>& x,
          const distance<TRep2, TInterval2>& y)
{
    using common_type = typename SIDIM_STD::common_type<
            distance<TRep1, TInterval1>,
            distance<TRep2, TInterval2> >::type;
    return common_type(common_type(x).count() - common_type(y).count());
}

// TODO: operator*
// TODO: operator/
// TODO: operator%

// ----=====================================================================----
//     distance_cast
// ----=====================================================================----

namespace sidim_detail
{

template <typename TFromDistance, typename TToDistance, typename RatioT,
          bool RatioNumeratorEqualsOne, bool RatioDenominatorEqualsOne>
struct distance_cast_helper;

// Special case R = 1.
template <typename TFromDistance, typename TToDistance, typename RatioT>
struct distance_cast_helper<TFromDistance, TToDistance, RatioT, true, true>
{
    constexpr TToDistance cast(const TFromDistance& from) const
    {
        return TToDistance(static_cast<typename TToDistance::rep>(
                               from.count()));
    }
};

// Special case R = rN / 1, rN != 1.
template <typename TFromDistance, typename TToDistance, typename RatioT>
struct distance_cast_helper<TFromDistance, TToDistance, RatioT, false, true>
{
    constexpr TToDistance cast(const TFromDistance& from) const
    {
        using common_type = typename SIDIM_STD::common_type<
                    typename TFromDistance::rep,
                    typename TToDistance::rep,
                    intmax_t>::type;

        return TToDistance(static_cast<typename TToDistance::rep>(
                               static_cast<common_type>(from.count())
                               * static_cast<common_type>(RatioT::num)));
    }
};

// Special case R = 1 / rD, rD != 1.
template <typename TFromDistance, typename TToDistance, typename RatioT>
struct distance_cast_helper<TFromDistance, TToDistance, RatioT, true, false>
{
    constexpr TToDistance cast(const TFromDistance& from) const
    {
        using common_type = typename SIDIM_STD::common_type<
                    typename TFromDistance::rep,
                    typename TToDistance::rep,
                    intmax_t>::type;

        return TToDistance(static_cast<typename TToDistance::rep>(
                               static_cast<common_type>(from.count())
                               / static_cast<common_type>(RatioT::den)));
    }
};

// General case R = rN / rD, rN != 1, rD != 1.
template <typename TFromDistance, typename TToDistance, typename RatioT>
struct distance_cast_helper<TFromDistance, TToDistance, RatioT, false, false>
{
    constexpr TToDistance cast(const TFromDistance& from) const
    {
        using common_type = typename SIDIM_STD::common_type<
                    typename TFromDistance::rep,
                    typename TToDistance::rep,
                    intmax_t>::type;

        return TToDistance(static_cast<typename TToDistance::rep>(
                               static_cast<common_type>(from.count())
                               * static_cast<common_type>(RatioT::num)
                               / static_cast<common_type>(RatioT::den)));
    }
};

template <typename TFromDistance, typename TToDistance>
struct distance_caster
{
    typedef typename SIDIM_STD::ratio_divide<
                         typename TFromDistance::interval,
                         typename TToDistance::interval>::type ratio;

    typedef distance_cast_helper<TFromDistance, TToDistance,
                                 ratio,
                                 ratio::num == 1,
                                 ratio::den == 1> helper;

    constexpr
    TToDistance cast(const TFromDistance& from) const
    {
        return helper().cast(from);
    }
};

} // namespace sidim_detail

//! A utility function to cast distances.
//! Cast from a <tt>distance<TRep, TInterval></tt> given in \p d to another
//! distance templated by \p TToDistance. The call
//! <tt>distance_cast<T>(d)</tt> is equivalent to
//! <tt>d.count() * d::interval / T::interval</tt>. If the destination interval
//! is coarser than the source interval, a truncation occurs if the destination
//! representation is not a floating point type.
//! All values are cast to at least intmax_t before performing the computation.
//!
//! The following example demonstrates, how to cast a distance in micrometers
//! to millimeters and meters:
//! \code
//! micrometer um(1234567);
//! auto mm = distance_cast<millimeter>(um); // mm.count() == 1234
//! auto m = distance_cast<meter>(m);        // m.count() == 1
//! \endcode
template <typename TToDistance, typename TRep, typename TInterval>
inline constexpr
typename SIDIM_STD::enable_if<sidim_detail::is_distance<TToDistance>::value,
                             TToDistance>::type
distance_cast(const distance<TRep, TInterval>& d)
{
    return sidim_detail::distance_caster<distance<TRep, TInterval>,
                                        TToDistance>().cast(d);
}

// ----=====================================================================----
//     position
// ----=====================================================================----

//! A position.
template <typename TTag, typename TDistance>
class position
{
    static_assert(sidim_detail::is_distance<TDistance>::value,
                  "The second template parameter must be a distance");

public:
    typedef TTag tag;
    typedef TDistance distance;
    typedef typename distance::rep rep;
    typedef typename distance::interval interval;

    constexpr
    position()
        : m_distance(distance::zero())
    {
    }

    //! Creates a position from a distance.
    //! Creates a position which at a distance \p d from the origin.
    constexpr explicit
    position(const distance& d)
        : m_distance(d)
    {
    }

    template <typename TDistance2,
              typename = typename SIDIM_STD::enable_if<
                  SIDIM_STD::is_convertible<TDistance2, distance>::value>::type>
    constexpr
    position(const position<tag, TDistance2>& pos)
        : m_distance(pos.distance_from_origin())
    {
    }

    //! Returns the distance from the origin.
    constexpr
    distance distance_from_origin() const
    {
        return m_distance;
    }

    // Arithmetic operators.

    //! Adds a distance.
    position& operator+=(const distance& d)
    {
        m_distance += d;
        return *this;
    }

    //! Subtracts a distance.
    position& operator-=(const distance& d)
    {
        m_distance -= d;
        return *this;
    }

    // Special values.

    static constexpr
    position max()
    {
        return position(distance::max());
    }

    static constexpr
    position min()
    {
        return position(distance::min());
    }

private:
    distance m_distance;
};

} // namespace sidim

// ----=====================================================================----
//     Specialization of common_type for sidim::position<>
// ----=====================================================================----

namespace SIDIM_STD
{

template <typename TTag, typename TDistance1, typename TDistance2>
struct common_type<sidim::position<TTag, TDistance1>,
                   sidim::position<TTag, TDistance2> >
{
    using type =
        sidim::position<TTag, typename common_type<TDistance1, TDistance2>::type>;
};

} // namespace SIDIM_STD

namespace sidim
{

// ----=====================================================================----
//     position_cast
// ----=====================================================================----

template <typename TToDistance, typename TTag, typename TFromDistance>
inline constexpr
position<TTag, TToDistance>
position_cast(const position<TTag, TFromDistance>& pos)
{
    return position<TTag, TToDistance>(distance_cast<TToDistance>(
                            pos.distance_from_origin()));
}

// ----=====================================================================----
//     position comparisons
// ----=====================================================================----

template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
bool
operator==(const position<TTag, TDistance1>& x,
           const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() == y.distance_from_origin();
}

template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
bool
operator!=(const position<TTag, TDistance1>& x,
           const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() != y.distance_from_origin();
}

template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
bool
operator<(const position<TTag, TDistance1>& x,
          const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() < y.distance_from_origin();
}

template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
bool
operator>(const position<TTag, TDistance1>& x,
          const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() > y.distance_from_origin();
}


template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
bool
operator<=(const position<TTag, TDistance1>& x,
           const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() <= y.distance_from_origin();
}

template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
bool
operator>=(const position<TTag, TDistance1>& x,
           const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() >= y.distance_from_origin();
}

// ----=====================================================================----
//     position arithmetic
// ----=====================================================================----

template <typename TTag, typename TDistance1,
          typename TRep2, typename TInterval2>
inline constexpr
position<TTag,
         typename SIDIM_STD::common_type<TDistance1, distance<TRep2, TInterval2> >::type>
operator+(const position<TTag, TDistance1>& pos,
          const distance<TRep2, TInterval2>& d)
{
    using common_distance = typename SIDIM_STD::common_type<
                                TDistance1, distance<TRep2, TInterval2> >::type;

    return position<TTag, common_distance>(pos.distance_from_origin() + d);
}

template <typename TRep1, typename TInterval1, typename TTag, typename TDistance2>
inline constexpr
position<TTag,
         typename SIDIM_STD::common_type<distance<TRep1, TInterval1>, TDistance2>::type>
operator+(const distance<TRep1, TInterval1>& d,
          const position<TTag, TDistance2>& pos)
{
    return pos + d;
}

template <typename TTag, typename TDistance1, typename TRep2, typename TInterval2>
inline constexpr
position<TTag,
         typename SIDIM_STD::common_type<TDistance1, distance<TRep2, TInterval2> >::type>
operator-(const position<TTag, TDistance1>& pos,
          const distance<TRep2, TInterval2>& d)
{
    using common_distance = typename SIDIM_STD::common_type<
                                TDistance1, distance<TRep2, TInterval2> >::type;

    return position<TTag, common_distance>(pos.distance_from_origin() - d);
}

template <typename TTag, typename TDistance1, typename TDistance2>
inline constexpr
typename SIDIM_STD::common_type<TDistance1, TDistance2>::type
operator-(const position<TTag, TDistance1>& x,
          const position<TTag, TDistance2>& y)
{
    return x.distance_from_origin() - y.distance_from_origin();
}

} // namespace sidim

#endif // SIDIM_HPP
