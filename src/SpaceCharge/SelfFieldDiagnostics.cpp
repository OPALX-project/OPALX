/**
 * @file SelfFieldDiagnostics.cpp
 * @brief Implements self-field lifecycle event recording.
 */

#include "SpaceCharge/SelfFieldDiagnostics.h"

#include "Utilities/OpalException.h"

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

    SelfFieldDiagnostics::SelfFieldDiagnostics(
            SelfFieldDiagnosticSchedule schedule, SelfFieldDiagnosticSink* sink)
        : sink_m(sink), schedule_m(schedule) {}

    SelfFieldDiagnostics::ScopedEvent SelfFieldDiagnostics::scopedEvent(
            SelfFieldEventKind kind, std::string_view label) {
        if (kind == SelfFieldEventKind::Count) {
            throw OpalException(
                    "SelfFieldDiagnostics::scopedEvent",
                    "SelfFieldEventKind::Count is not a lifecycle event.");
        }
        return ScopedEvent(*this, kind, label);
    }

    std::size_t SelfFieldDiagnostics::completedEventCount(SelfFieldEventKind kind) const {
        if (kind == SelfFieldEventKind::Count) {
            throw OpalException(
                    "SelfFieldDiagnostics::completedEventCount",
                    "SelfFieldEventKind::Count is not a lifecycle event.");
        }
        return completedEvents_m[index(kind)];
    }

    std::chrono::nanoseconds SelfFieldDiagnostics::totalEventDuration(
            SelfFieldEventKind kind) const {
        if (kind == SelfFieldEventKind::Count) {
            throw OpalException(
                    "SelfFieldDiagnostics::totalEventDuration",
                    "SelfFieldEventKind::Count is not a lifecycle event.");
        }
        return totalDurations_m[index(kind)];
    }

    bool SelfFieldDiagnostics::shouldPrintBinTable(long long step) const noexcept {
        return step >= 0 && schedule_m.binTableFrequency > 0
               && static_cast<std::size_t>(step) % schedule_m.binTableFrequency == 0;
    }

    bool SelfFieldDiagnostics::shouldDumpPlane(long long step) const noexcept {
        return step >= 0 && schedule_m.planeDumpFrequency > 0
               && static_cast<std::size_t>(step) % schedule_m.planeDumpFrequency == 0;
    }

    void SelfFieldDiagnostics::reset() {
        completedEvents_m.fill(0);
        totalDurations_m.fill(std::chrono::nanoseconds::zero());
        backendSolveCount_m   = 0;
        redistributionCount_m = 0;
        currentBinCount_m     = 0;
        hasCurrentBinCount_m  = false;
    }

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
        totalDurations_m[index(kind)] += duration;
        if (sink_m != nullptr) {
            sink_m->eventFinished(kind, label, duration);
        }
    }

}  // namespace opalx::spacecharge
