#ifndef TPS_HOST_STORAGE_HPP
#define TPS_HOST_STORAGE_HPP
#include <Kokkos_Core.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
namespace tps::detail {
    /** @brief Kokkos-owned host allocation with host-managed C++ element lifetimes.
     * @details WithoutInitializing prevents nested View-owning objects from being
     * constructed/destructed inside Kokkos allocation callbacks. All elements are
     * value-initialized on the host and destroyed before the allocation is released.
     */
    template <class T>
    class HostStorage {
    public:
        HostStorage() = default;
        HostStorage(const char* label, int size) : memory_(allocate(label, size)), size_(size) {
            std::uninitialized_value_construct_n(memory_.data(), size_);
        }
        HostStorage(const HostStorage&)            = delete;
        HostStorage& operator=(const HostStorage&) = delete;
        HostStorage(HostStorage&& other) noexcept
            : memory_(std::move(other.memory_)), size_(std::exchange(other.size_, 0)) {}
        HostStorage& operator=(HostStorage&& other) noexcept {
            if (this != &other) {
                clear();
                memory_ = std::move(other.memory_);
                size_   = std::exchange(other.size_, 0);
            }
            return *this;
        }
        ~HostStorage() { clear(); }
        T* data() const noexcept { return memory_.data(); }

    private:
        static auto allocate(const char* label, int size) {
            if (!Kokkos::is_initialized() || Kokkos::is_finalized())
                throw std::logic_error("TPS arrays require an active Kokkos runtime");
            return Kokkos::View<T*, Kokkos::HostSpace>(
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(label)), size);
        }
        void clear() noexcept {
            if (size_) std::destroy_n(memory_.data(), std::exchange(size_, 0));
        }
        Kokkos::View<T*, Kokkos::HostSpace> memory_;
        int size_ = 0;
    };
}  // namespace tps::detail
#endif
