/**
 * @file SelfFieldDiagnostics.cpp
 * @brief Implements self-field lifecycle event recording.
 */

#include "SpaceCharge/SelfFieldDiagnostics.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace opalx::spacecharge {

    SelfFieldDiagnostics::ScopedEvent::ScopedEvent(
            SelfFieldDiagnostics& diagnostics, SelfFieldEventKind kind, std::string_view label)
        : diagnostics_m(&diagnostics), kind_m(kind), label_m(label), start_m(clock_type::now()) {
        diagnostics_m->eventStarted(kind_m, label_m);
    }

    SelfFieldDiagnostics::ScopedEvent::~ScopedEvent() { finish(); }

    SelfFieldDiagnostics::ScopedEvent::ScopedEvent(ScopedEvent&& other) noexcept
        : diagnostics_m(std::exchange(other.diagnostics_m, nullptr)),
          kind_m(other.kind_m),
          label_m(other.label_m),
          start_m(other.start_m) {}

    SelfFieldDiagnostics::ScopedEvent& SelfFieldDiagnostics::ScopedEvent::operator=(
            ScopedEvent&& other) noexcept {
        if (this != &other) {
            finish();
            diagnostics_m = std::exchange(other.diagnostics_m, nullptr);
            kind_m        = other.kind_m;
            label_m       = other.label_m;
            start_m       = other.start_m;
        }
        return *this;
    }

    void SelfFieldDiagnostics::ScopedEvent::finish() noexcept {
        if (diagnostics_m == nullptr) {
            return;
        }
        const auto duration =
                std::chrono::duration_cast<std::chrono::nanoseconds>(clock_type::now() - start_m);
        diagnostics_m->eventFinished(kind_m, label_m, duration);
        diagnostics_m = nullptr;
    }

    SelfFieldDiagnostics::SelfFieldDiagnostics(SelfFieldDiagnosticSink* sink) : sink_m(sink) {}

    SelfFieldDiagnostics::ScopedEvent SelfFieldDiagnostics::scopedEvent(
            SelfFieldEventKind kind, std::string_view label) {
        if (kind == SelfFieldEventKind::Count) {
            throw std::invalid_argument("SelfFieldEventKind::Count is not a lifecycle event");
        }
        return ScopedEvent(*this, kind, label);
    }

    std::size_t SelfFieldDiagnostics::completedEventCount(SelfFieldEventKind kind) const {
        if (kind == SelfFieldEventKind::Count) {
            throw std::invalid_argument("SelfFieldEventKind::Count is not a lifecycle event");
        }
        return completedEvents_m[index(kind)];
    }

    void SelfFieldDiagnostics::reset() { completedEvents_m.fill(0); }

    void SelfFieldDiagnostics::eventStarted(
            SelfFieldEventKind kind, std::string_view label) noexcept {
        if (sink_m != nullptr) {
            sink_m->eventStarted(kind, label);
        }
    }

    void SelfFieldDiagnostics::eventFinished(
            SelfFieldEventKind kind, std::string_view label,
            std::chrono::nanoseconds duration) noexcept {
        ++completedEvents_m[index(kind)];
        if (sink_m != nullptr) {
            sink_m->eventFinished(kind, label, duration);
        }
    }

}  // namespace opalx::spacecharge
