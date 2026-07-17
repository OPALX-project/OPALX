/**
 * @file SelfFieldDiagnostics.h
 * @brief Stable host-side lifecycle diagnostics for self-field algorithms.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_DIAGNOSTICS_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_DIAGNOSTICS_H

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace opalx::spacecharge {

    /** @brief Stable lifecycle boundaries visible to self-field diagnostics. */
    enum class SelfFieldEventKind : std::uint8_t {
        Solve,
        DomainUpdate,
        SolveUnit,
        SolvePass,
        BackendSolve,
        FieldComposition,
        Count
    };

    /**
     * @brief Optional observer for lifecycle events.
     *
     * A sink records work chosen by the algorithm. It must not alter the solve plan or throw from
     * callbacks. Kernel-level timing remains outside this stable interface.
     */
    class SelfFieldDiagnosticSink {
    public:
        virtual ~SelfFieldDiagnosticSink() = default;

        virtual void eventStarted(SelfFieldEventKind kind, std::string_view label) noexcept = 0;
        virtual void eventFinished(
                SelfFieldEventKind kind, std::string_view label,
                std::chrono::nanoseconds duration) noexcept = 0;
    };

    /** @brief Lifecycle recorder passed to one concrete SelfFieldAlgorithm. */
    class SelfFieldDiagnostics {
    public:
        using clock_type = std::chrono::steady_clock;

        /** @brief Exception-safe scope for one lifecycle event. */
        class ScopedEvent {
        public:
            ScopedEvent(
                    SelfFieldDiagnostics& diagnostics, SelfFieldEventKind kind,
                    std::string_view label);
            ~ScopedEvent();

            ScopedEvent(const ScopedEvent&)            = delete;
            ScopedEvent& operator=(const ScopedEvent&) = delete;
            ScopedEvent(ScopedEvent&& other) noexcept;
            ScopedEvent& operator=(ScopedEvent&& other) noexcept;

            /** @brief Finish the event before the scope ends. Safe to call more than once. */
            void finish() noexcept;

        private:
            SelfFieldDiagnostics* diagnostics_m = nullptr;
            SelfFieldEventKind kind_m           = SelfFieldEventKind::Solve;
            std::string_view label_m;
            clock_type::time_point start_m;
        };

        explicit SelfFieldDiagnostics(SelfFieldDiagnosticSink* sink = nullptr);

        /** @brief Start a lifecycle event and return its exception-safe scope. */
        [[nodiscard]] ScopedEvent scopedEvent(SelfFieldEventKind kind, std::string_view label = {});

        /** @brief Number of successfully closed events of one kind. */
        [[nodiscard]] std::size_t completedEventCount(SelfFieldEventKind kind) const;

        /** @brief Reset lifecycle counters, for example after construction warm-up solves. */
        void reset();

        /** @brief Replace the optional non-owning sink. */
        void setSink(SelfFieldDiagnosticSink* sink) { sink_m = sink; }

    private:
        friend class ScopedEvent;

        static constexpr std::size_t index(SelfFieldEventKind kind) {
            return static_cast<std::size_t>(kind);
        }

        void eventStarted(SelfFieldEventKind kind, std::string_view label) noexcept;
        void eventFinished(
                SelfFieldEventKind kind, std::string_view label,
                std::chrono::nanoseconds duration) noexcept;

        SelfFieldDiagnosticSink* sink_m = nullptr;
        std::array<std::size_t, static_cast<std::size_t>(SelfFieldEventKind::Count)>
                completedEvents_m{};
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_DIAGNOSTICS_H
