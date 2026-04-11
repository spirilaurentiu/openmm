#ifndef OPENMM_CUDAARRAY_H_
#define OPENMM_CUDAARRAY_H_

#include "openmm/common/windowsExportCommon.h"
#include "openmm/common/ArrayInterface.h"
#include <cuda.h>
#include <memory> 
#include <string>
#include <vector>

namespace OpenMM {

class CudaContext;

/**
 * This class encapsulates a block of CUDA device memory. 
 * Refactored for C++17: Move-only RAII container.
 */
class OPENMM_EXPORT_COMMON CudaArray : public ArrayInterface {
public:
    /**
     * Factory method: Returns a unique_ptr to ensure the caller assumes ownership.
     */
    template <class T>
    [[nodiscard]] static std::unique_ptr<CudaArray> create(CudaContext& context, size_t size, const std::string& name) {
        return std::make_unique<CudaArray>(context, size, sizeof(T), name);
    }

    CudaArray();
    CudaArray(CudaContext& context, size_t size, int elementSize, const std::string& name);
    
    /**
     * RAII Destructor: Marked noexcept to ensure safety during stack unwinding.
     */
    ~CudaArray() noexcept;

    // --- Rule of Five ---
    // CUDA memory is a unique resource; we allow moving, but forbid copying.
    CudaArray(const CudaArray&) = delete;
    CudaArray& operator=(const CudaArray&) = delete;
    CudaArray(CudaArray&& other) noexcept;
    CudaArray& operator=(CudaArray&& other) noexcept;

    void initialize(ComputeContext& context, size_t size, int elementSize, const std::string& name);

    template <class T>
    void initialize(ComputeContext& context, size_t size, const std::string& name) {
        initialize(context, size, sizeof(T), name);
    }

    void resize(size_t size);

    [[nodiscard]] bool isInitialized() const noexcept {
        return (pointer != 0);
    }

    [[nodiscard]] size_t getSize() const noexcept {
        return size;
    }

    [[nodiscard]] int getElementSize() const noexcept {
        return elementSize;
    }

    [[nodiscard]] const std::string& getName() const noexcept {
        return name;
    }

    [[nodiscard]] ComputeContext& getContext();

    [[nodiscard]] CUdeviceptr& getDevicePointer() noexcept {
        return pointer;
    }

    // --- Data Transfer ---

    template <class T>
    void upload(const std::vector<T>& data, bool convert = false) {
        ArrayInterface::upload(data, convert);
    }

    template <class T>
    void download(std::vector<T>& data) const {
        ArrayInterface::download(data);
    }

    void upload(const void* data, bool blocking = true) {
        uploadSubArray(data, 0, static_cast<int>(getSize()), blocking);
    }

    void uploadSubArray(const void* data, int offset, int elements, bool blocking = true);
    
    void download(void* data, bool blocking = true) const;
    
    void copyTo(ArrayInterface& dest) const;

private:
    CudaContext* context {nullptr};
    CUdeviceptr pointer {0};
    size_t size {0};
    int elementSize {0};
    bool ownsMemory {false};
    std::string name;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAARRAY_H_*/