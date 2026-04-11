#include "CudaArray.h"
#include "CudaContext.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/OpenMMException.h"
#include <iostream>
#include <string_view>
#include <utility>

using namespace OpenMM;

// Helper to handle CUDA errors without repetitive boilerplate
void checkCuda(CUresult result, std::string_view action, const std::string& name) {
    if (result != CUDA_SUCCESS) {
        throw OpenMMException(std::string("Error ") + action.data() + " array " + name + 
                             ": " + CudaContext::getErrorString(result) + 
                             " (" + std::to_string(result) + ")");
    }
}

CudaArray::CudaArray() : pointer(0), context(nullptr), size(0), elementSize(0), ownsMemory(false) {}

CudaArray::CudaArray(CudaContext& context, size_t size, int elementSize, const std::string& name) 
    : pointer(0), context(nullptr), ownsMemory(false) {
    initialize(context, size, elementSize, name);
}

// Move Constructor
CudaArray::CudaArray(CudaArray&& other) noexcept 
    : pointer(other.pointer), context(other.context), size(other.size), 
      elementSize(other.elementSize), name(std::move(other.name)), ownsMemory(other.ownsMemory) {
    other.pointer = 0;
    other.ownsMemory = false;
}

// Move Assignment
CudaArray& CudaArray::operator=(CudaArray&& other) noexcept {
    if (this != &other) {
        // Clean up existing resource
        if (pointer != 0 && ownsMemory && context && context->getContextIsValid()) {
            ContextSelector selector(*context);
            cuMemFree(pointer);
        }
        
        pointer = other.pointer;
        context = other.context;
        size = other.size;
        elementSize = other.elementSize;
        name = std::move(other.name);
        ownsMemory = other.ownsMemory;

        other.pointer = 0;
        other.ownsMemory = false;
    }
    return *this;
}

CudaArray::~CudaArray() noexcept {
    if (pointer != 0 && ownsMemory && context && context->getContextIsValid()) {
        try {
            ContextSelector selector(*context);
            CUresult result = cuMemFree(pointer);
            if (result != CUDA_SUCCESS) {
                // Destructors must not throw. Log the error instead.
                std::cerr << "Error deleting array " << name << ": " 
                          << CudaContext::getErrorString(result) << std::endl;
            }
        } catch (...) {
            // Catch all to prevent exceptions escaping the destructor
        }
    }
}

void CudaArray::initialize(ComputeContext& context, size_t size, int elementSize, const std::string& name) {
    if (this->pointer != 0) {
        throw OpenMMException("CudaArray has already been initialized");
    }

    this->context = &dynamic_cast<CudaContext&>(context);
    this->size = size;
    this->elementSize = elementSize;
    this->name = name;
    this->ownsMemory = true;

    ContextSelector selector(*this->context);
    checkCuda(cuMemAlloc(&pointer, size * elementSize), "creating", name);
}

void CudaArray::resize(size_t newSize) {
    if (pointer == 0) {
        throw OpenMMException("CudaArray has not been initialized");
    }
    if (!ownsMemory) {
        throw OpenMMException("Cannot resize an array that does not own its storage");
    }

    {
        ContextSelector selector(*context);
        checkCuda(cuMemFree(pointer), "deleting", name);
    }

    pointer = 0;
    initialize(*context, newSize, elementSize, name);
}

[[nodiscard]] ComputeContext& CudaArray::getContext() {
    return *context;
}

void CudaArray::uploadSubArray(const void* data, int offset, int elements, bool blocking) {
    if (pointer == 0) {
        throw OpenMMException("CudaArray has not been initialized");
    }
    if (offset < 0 || offset + elements > static_cast<int>(size)) {
        throw OpenMMException("uploadSubArray: data exceeds range of array");
    }

    ContextSelector selector(*context);
    CUresult result;
    size_t byteOffset = static_cast<size_t>(offset) * elementSize;
    size_t byteCount = static_cast<size_t>(elements) * elementSize;

    if (blocking) {
        result = cuMemcpyHtoD(pointer + byteOffset, data, byteCount);
    } else {
        result = cuMemcpyHtoDAsync(pointer + byteOffset, data, byteCount, context->getCurrentStream());
    }
    
    checkCuda(result, "uploading", name);
}

void CudaArray::download(void* data, bool blocking) const {
    if (pointer == 0) {
        throw OpenMMException("CudaArray has not been initialized");
    }

    ContextSelector selector(*context);
    CUresult result;
    if (blocking) {
        result = cuMemcpyDtoH(data, pointer, size * elementSize);
    } else {
        result = cuMemcpyDtoHAsync(data, pointer, size * elementSize, context->getCurrentStream());
    }

    checkCuda(result, "downloading", name);
}

void CudaArray::copyTo(ArrayInterface& dest) const {
    if (pointer == 0) {
        throw OpenMMException("CudaArray has not been initialized");
    }
    if (dest.getSize() != size || dest.getElementSize() != elementSize) {
        throw OpenMMException("Error copying array " + name + " to " + dest.getName() + 
                             ": The destination array does not match the size of the array");
    }

    CudaArray& cuDest = context->unwrap(dest);
    ContextSelector selector(*context);
    checkCuda(cuMemcpyDtoDAsync(cuDest.getDevicePointer(), pointer, size * elementSize, context->getCurrentStream()),
              "copying", name);
}