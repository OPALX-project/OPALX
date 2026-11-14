/**
 * @file MomentRequest.h
 * @brief Groups of beam moments that a caller can request.
 */

#ifndef OPALX_MOMENTS_MOMENT_REQUEST_H
#define OPALX_MOMENTS_MOMENT_REQUEST_H

#include <cstdint>

namespace opalx::moments {

    /**
     * @brief One group of beam moments, the unit of a request.
     *
     * Groups follow the particle passes: Bounds, MeanPosition and MeanMomentum share one
     * first-order pass that reads only the arrays they need, Spread is a second pass centred on
     * both means, and Polarization has its own pass.
     */
    enum class Moment : std::uint8_t {
        Bounds       = 1U << 0U,  ///< Smallest and largest position per axis.
        MeanPosition = 1U << 1U,  ///< Mean position.
        MeanMomentum = 1U << 2U,  ///< Mean beta*gamma, mean gamma and mean kinetic energy.
        /// Central 6x6 covariance, kinetic-energy spread and rest-frame temperature.
        Spread       = 1U << 3U,
        Polarization = 1U << 4U,  ///< Mean and rms polarization; only for containers with spin.
    };

    /** @brief Set of moment groups; a single Moment converts to it implicitly. */
    class Request {
    public:
        constexpr Request() = default;
        constexpr Request(Moment moment) : bits_m(static_cast<std::uint8_t>(moment)) {}

        [[nodiscard]] constexpr bool contains(Request other) const {
            return (bits_m & other.bits_m) == other.bits_m;
        }
        [[nodiscard]] constexpr bool empty() const { return bits_m == 0; }
        [[nodiscard]] constexpr Request without(Request other) const {
            return fromBits(bits_m & ~other.bits_m);
        }

        friend constexpr Request operator|(Request left, Request right) {
            return fromBits(left.bits_m | right.bits_m);
        }
        constexpr Request& operator|=(Request other) { return *this = *this | other; }
        constexpr bool operator==(const Request&) const = default;

    private:
        static constexpr Request fromBits(unsigned bits) {
            Request request;
            request.bits_m = static_cast<std::uint8_t>(bits);
            return request;
        }

        std::uint8_t bits_m = 0;
    };

    constexpr Request operator|(Moment left, Moment right) { return Request(left) | right; }

    /// Every group, as read by the statistics output.
    inline constexpr Request allMoments = Moment::Bounds | Moment::MeanPosition
                                          | Moment::MeanMomentum | Moment::Spread
                                          | Moment::Polarization;

}  // namespace opalx::moments

#endif  // OPALX_MOMENTS_MOMENT_REQUEST_H
