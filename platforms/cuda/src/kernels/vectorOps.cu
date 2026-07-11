/**
 * @file vectorOps.cu
 * @brief Component-wise arithmetic and geometric helpers for the built-in CUDA
 *        vector types (int/float/double, widths 2/3/4; plus short4 in trimTo3).
 *
 * All symbols are pure inline __device__ helpers with no kernels: single-scalar
 * make_* constructors, operator overloads, and the named helpers dot, cross,
 * normalize, and trimTo3. They depend only on the CUDA runtime's base make_TN
 * constructors and rsqrt/rsqrtf; nothing here is injected by NVRTC macros.
 *
 * @par Execution context
 * Every function is a leaf __device__ helper: no shared memory, no
 * synchronization, not warp-synchronous, no global or static state. Each thread
 * evaluates its call independently, so all are safe under arbitrary warp
 * divergence and impose no participation requirement. Results are deterministic
 * functions of the inputs (no atomics, no cross-thread communication).
 *
 * @par Numerical contract (holds file-wide unless a function states otherwise)
 * - Semantics are COMPONENT-WISE (Hadamard): operator* is never a dot or cross
 *   product, and there is no matrix algebra.
 * - Overloads are element-type- and width-HOMOGENEOUS: both operands and the
 *   result share one element type and one width. No mixed int/float/double or
 *   mixed-width overload exists, so implicit scalar promotion never engages; a
 *   mixed-type call fails to compile rather than silently converting.
 * - No divide-by-zero or zero-length guarding is performed anywhere: integer
 *   division by zero is undefined behavior; float/double division and normalize
 *   of a zero-length vector yield inf/nan rather than trapping.
 * - Operands are passed by value and results returned by value (compound
 *   assignment excepted); there are no pointer parameters and no aliasing
 *   hazard, and passing one object as both operands is well defined.
 * - Products written as (a*b)+(c*d) are subject to the compiler's FMA
 *   contraction (nvcc --fmad): contraction fuses a multiply into an adjacent add
 *   and can change the last ulp, so the exact bit patterns of dot, cross, and
 *   normalize depend on the JIT's fmad setting.
 */

/**
 * @defgroup make_broadcast Single-argument make_* broadcast constructors
 * @brief Build a vector by replicating one scalar into every component.
 *
 * Each overload forwards to the built-in multi-argument constructor, so every
 * component of the result equals @p value. The argument's type selects the
 * overload and fixes the result element type; no int/float conversion occurs.
 * Covers int/float/double in widths 2/3/4.
 *
 * @param[in] value scalar copied unchanged into every component.
 * @return vector whose components all equal @p value.
 * @{
 */
// Versions of make_x() that take a single value and set all components to that.

inline __device__ auto make_int2(int value) -> int2 {
    return make_int2(value, value);
}

inline __device__ auto make_int3(int value) -> int3 {
    return make_int3(value, value, value);
}

inline __device__ auto make_int4(int value) -> int4 {
    return make_int4(value, value, value, value);
}

inline __device__ auto make_float2(float value) -> float2 {
    return make_float2(value, value);
}

inline __device__ auto make_float3(float value) -> float3 {
    return make_float3(value, value, value);
}

inline __device__ auto make_float4(float value) -> float4 {
    return make_float4(value, value, value, value);
}

inline __device__ auto make_double2(double value) -> double2 {
    return make_double2(value, value);
}

inline __device__ auto make_double3(double value) -> double3 {
    return make_double3(value, value, value);
}

inline __device__ auto make_double4(double value) -> double4 {
    return make_double4(value, value, value, value);
}

/** @} */

/**
 * @defgroup vec_negate Component-wise unary negation
 * @brief Negate every component independently: result.c = -value.c.
 *
 * The unary operator- for int/float/double, widths 2/3/4. Arity distinguishes it
 * from the binary vector-vector subtraction below.
 *
 * @param[in] value operand.
 * @return the component-wise negation of @p value.
 * @note For the int overloads, negating the minimum representable value is
 *       signed-integer overflow (undefined behavior); float/double negation is
 *       exact and never overflows.
 * @{
 */
// Negate a vector.

inline __device__ auto operator-(int2 value) -> int2 {
    return make_int2(-value.x, -value.y);
}

inline __device__ auto operator-(int3 value) -> int3 {
    return make_int3(-value.x, -value.y, -value.z);
}

inline __device__ auto operator-(int4 value) -> int4 {
    return make_int4(-value.x, -value.y, -value.z, -value.w);
}

inline __device__ auto operator-(float2 value) -> float2 {
    return make_float2(-value.x, -value.y);
}

inline __device__ auto operator-(float3 value) -> float3 {
    return make_float3(-value.x, -value.y, -value.z);
}

inline __device__ auto operator-(float4 value) -> float4 {
    return make_float4(-value.x, -value.y, -value.z, -value.w);
}

inline __device__ auto operator-(double2 value) -> double2 {
    return make_double2(-value.x, -value.y);
}

inline __device__ auto operator-(double3 value) -> double3 {
    return make_double3(-value.x, -value.y, -value.z);
}

inline __device__ auto operator-(double4 value) -> double4 {
    return make_double4(-value.x, -value.y, -value.z, -value.w);
}

/** @} */

/**
 * @defgroup vec_elementwise Component-wise binary arithmetic (vector op vector)
 * @brief Hadamard +, -, *, / of two same-type vectors: result.c = lhs.c op rhs.c.
 *
 * operator* is the element-wise product, not a dot or cross product. Covers
 * int/float/double, widths 2/3/4, for all four operators.
 *
 * @param[in] lhs left operand.
 * @param[in] rhs right operand.
 * @return the component-wise result.
 * @note Division is per-component with no reciprocal shortcut: the float/double
 *       overloads perform a true IEEE divide on each component, in contrast to
 *       the vector/scalar division family, which multiplies by a precomputed
 *       reciprocal. Integer division truncates toward zero.
 * @warning A zero divisor is undefined behavior for the int overloads and yields
 *          inf/nan for float/double; no divisor check is performed.
 * @{
 */
// Add two vectors.

inline __device__ auto operator+(int2 lhs, int2 rhs) -> int2 {
    return make_int2(lhs.x+rhs.x, lhs.y+rhs.y);
}

inline __device__ auto operator+(int3 lhs, int3 rhs) -> int3 {
    return make_int3(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
}

inline __device__ auto operator+(int4 lhs, int4 rhs) -> int4 {
    return make_int4(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z, lhs.w+rhs.w);
}

inline __device__ auto operator+(float2 lhs, float2 rhs) -> float2 {
    return make_float2(lhs.x+rhs.x, lhs.y+rhs.y);
}

inline __device__ auto operator+(float3 lhs, float3 rhs) -> float3 {
    return make_float3(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
}

inline __device__ auto operator+(float4 lhs, float4 rhs) -> float4 {
    return make_float4(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z, lhs.w+rhs.w);
}

inline __device__ auto operator+(double2 lhs, double2 rhs) -> double2 {
    return make_double2(lhs.x+rhs.x, lhs.y+rhs.y);
}

inline __device__ auto operator+(double3 lhs, double3 rhs) -> double3 {
    return make_double3(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
}

inline __device__ auto operator+(double4 lhs, double4 rhs) -> double4 {
    return make_double4(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z, lhs.w+rhs.w);
}

// Subtract two vectors.

inline __device__ auto operator-(int2 lhs, int2 rhs) -> int2 {
    return make_int2(lhs.x-rhs.x, lhs.y-rhs.y);
}

inline __device__ auto operator-(int3 lhs, int3 rhs) -> int3 {
    return make_int3(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}

inline __device__ auto operator-(int4 lhs, int4 rhs) -> int4 {
    return make_int4(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z, lhs.w-rhs.w);
}

inline __device__ auto operator-(float2 lhs, float2 rhs) -> float2 {
    return make_float2(lhs.x-rhs.x, lhs.y-rhs.y);
}

inline __device__ auto operator-(float3 lhs, float3 rhs) -> float3 {
    return make_float3(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}

inline __device__ auto operator-(float4 lhs, float4 rhs) -> float4 {
    return make_float4(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z, lhs.w-rhs.w);
}

inline __device__ auto operator-(double2 lhs, double2 rhs) -> double2 {
    return make_double2(lhs.x-rhs.x, lhs.y-rhs.y);
}

inline __device__ auto operator-(double3 lhs, double3 rhs) -> double3 {
    return make_double3(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}

inline __device__ auto operator-(double4 lhs, double4 rhs) -> double4 {
    return make_double4(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z, lhs.w-rhs.w);
}

// Multiply two vectors.

inline __device__ auto operator*(int2 lhs, int2 rhs) -> int2 {
    return make_int2(lhs.x*rhs.x, lhs.y*rhs.y);
}

inline __device__ auto operator*(int3 lhs, int3 rhs) -> int3 {
    return make_int3(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
}

inline __device__ auto operator*(int4 lhs, int4 rhs) -> int4 {
    return make_int4(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z, lhs.w*rhs.w);
}

inline __device__ auto operator*(float2 lhs, float2 rhs) -> float2 {
    return make_float2(lhs.x*rhs.x, lhs.y*rhs.y);
}

inline __device__ auto operator*(float3 lhs, float3 rhs) -> float3 {
    return make_float3(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
}

inline __device__ auto operator*(float4 lhs, float4 rhs) -> float4 {
    return make_float4(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z, lhs.w*rhs.w);
}

inline __device__ auto operator*(double2 lhs, double2 rhs) -> double2 {
    return make_double2(lhs.x*rhs.x, lhs.y*rhs.y);
}

inline __device__ auto operator*(double3 lhs, double3 rhs) -> double3 {
    return make_double3(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
}

inline __device__ auto operator*(double4 lhs, double4 rhs) -> double4 {
    return make_double4(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z, lhs.w*rhs.w);
}

// Divide two vectors.

inline __device__ auto operator/(int2 lhs, int2 rhs) -> int2 {
    return make_int2(lhs.x/rhs.x, lhs.y/rhs.y);
}

inline __device__ auto operator/(int3 lhs, int3 rhs) -> int3 {
    return make_int3(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z);
}

inline __device__ auto operator/(int4 lhs, int4 rhs) -> int4 {
    return make_int4(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z, lhs.w/rhs.w);
}

inline __device__ auto operator/(float2 lhs, float2 rhs) -> float2 {
    return make_float2(lhs.x/rhs.x, lhs.y/rhs.y);
}

inline __device__ auto operator/(float3 lhs, float3 rhs) -> float3 {
    return make_float3(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z);
}

inline __device__ auto operator/(float4 lhs, float4 rhs) -> float4 {
    return make_float4(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z, lhs.w/rhs.w);
}

inline __device__ auto operator/(double2 lhs, double2 rhs) -> double2 {
    return make_double2(lhs.x/rhs.x, lhs.y/rhs.y);
}

inline __device__ auto operator/(double3 lhs, double3 rhs) -> double3 {
    return make_double3(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z);
}

inline __device__ auto operator/(double4 lhs, double4 rhs) -> double4 {
    return make_double4(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z, lhs.w/rhs.w);
}

/** @} */

/**
 * @defgroup vec_compound_vv Compound assignment (vector op= vector)
 * @brief In-place element-wise +=, -=, *=, /=: lhs.c op= rhs.c per component.
 *
 * @p lhs is taken by reference and mutated; @p rhs is taken by value. The return
 * type is void, so these do not chain (a += b += c does not compile). Covers
 * int/float/double, widths 2/3/4. Because @p rhs is copied, passing the same
 * object as both operands (a += a) is well defined.
 *
 * @param[in,out] lhs updated in place with the result.
 * @param[in]     rhs right operand.
 * @note Division semantics match the binary vector/vector family: integer
 *       truncation, and inf/nan (float/double) or undefined behavior (int) on a
 *       zero divisor. The scalar counterpart (vector op= scalar) is a separate
 *       family below.
 * @{
 */
// += operator

inline __device__ void operator+=(int2& lhs, int2 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y;
}

inline __device__ void operator+=(int3& lhs, int3 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z;
}

inline __device__ void operator+=(int4& lhs, int4 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; lhs.w += rhs.w;
}

inline __device__ void operator+=(float2& lhs, float2 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y;
}

inline __device__ void operator+=(float3& lhs, float3 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z;
}

inline __device__ void operator+=(float4& lhs, float4 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; lhs.w += rhs.w;
}

inline __device__ void operator+=(double2& lhs, double2 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y;
}

inline __device__ void operator+=(double3& lhs, double3 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z;
}

inline __device__ void operator+=(double4& lhs, double4 rhs) {
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; lhs.w += rhs.w;
}

// -= operator

inline __device__ void operator-=(int2& lhs, int2 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y;
}

inline __device__ void operator-=(int3& lhs, int3 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z;
}

inline __device__ void operator-=(int4& lhs, int4 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z; lhs.w -= rhs.w;
}

inline __device__ void operator-=(float2& lhs, float2 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y;
}

inline __device__ void operator-=(float3& lhs, float3 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z;
}

inline __device__ void operator-=(float4& lhs, float4 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z; lhs.w -= rhs.w;
}

inline __device__ void operator-=(double2& lhs, double2 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y;
}

inline __device__ void operator-=(double3& lhs, double3 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z;
}

inline __device__ void operator-=(double4& lhs, double4 rhs) {
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z; lhs.w -= rhs.w;
}

// *= operator

inline __device__ void operator*=(int2& lhs, int2 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y;
}

inline __device__ void operator*=(int3& lhs, int3 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z;
}

inline __device__ void operator*=(int4& lhs, int4 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z; lhs.w *= rhs.w;
}

inline __device__ void operator*=(float2& lhs, float2 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y;
}

inline __device__ void operator*=(float3& lhs, float3 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z;
}

inline __device__ void operator*=(float4& lhs, float4 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z; lhs.w *= rhs.w;
}

inline __device__ void operator*=(double2& lhs, double2 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y;
}

inline __device__ void operator*=(double3& lhs, double3 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z;
}

inline __device__ void operator*=(double4& lhs, double4 rhs) {
    lhs.x *= rhs.x; lhs.y *= rhs.y; lhs.z *= rhs.z; lhs.w *= rhs.w;
}

// /= operator

inline __device__ void operator/=(int2& lhs, int2 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y;
}

inline __device__ void operator/=(int3& lhs, int3 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z;
}

inline __device__ void operator/=(int4& lhs, int4 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z; lhs.w /= rhs.w;
}

inline __device__ void operator/=(float2& lhs, float2 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y;
}

inline __device__ void operator/=(float3& lhs, float3 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z;
}

inline __device__ void operator/=(float4& lhs, float4 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z; lhs.w /= rhs.w;
}

inline __device__ void operator/=(double2& lhs, double2 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y;
}

inline __device__ void operator/=(double3& lhs, double3 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z;
}

inline __device__ void operator/=(double4& lhs, double4 rhs) {
    lhs.x /= rhs.x; lhs.y /= rhs.y; lhs.z /= rhs.z; lhs.w /= rhs.w;
}

/** @} */

/**
 * @defgroup vec_scalar_mul Scalar-broadcast multiplication (vector * scalar)
 * @brief Scale every component by one scalar: result.c = vector.c * constant.
 *
 * Both operand orders, (vector, scalar) and (scalar, vector), are provided; they
 * are equal up to IEEE multiply commutativity. The scalar and vector element
 * types match per overload (int, float, double); there is no mixed-precision
 * form, so e.g. `2.0 * floatVec` does not resolve through these overloads.
 *
 * @param[in] vector   operand whose components are scaled.
 * @param[in] constant scalar multiplier.
 * @return the scaled vector.
 * @{
 */
// Multiply a vector by a constant.

inline __device__ auto operator*(int2 vector, int constant) -> int2 {
    return make_int2(vector.x*constant, vector.y*constant);
}

inline __device__ auto operator*(int3 vector, int constant) -> int3 {
    return make_int3(vector.x*constant, vector.y*constant, vector.z*constant);
}

inline __device__ auto operator*(int4 vector, int constant) -> int4 {
    return make_int4(vector.x*constant, vector.y*constant, vector.z*constant, vector.w*constant);
}

inline __device__ auto operator*(int constant, int2 vector) -> int2 {
    return make_int2(constant*vector.x, constant*vector.y);
}

inline __device__ auto operator*(int constant, int3 vector) -> int3 {
    return make_int3(constant*vector.x, constant*vector.y, constant*vector.z);
}

inline __device__ auto operator*(int constant, int4 vector) -> int4 {
    return make_int4(constant*vector.x, constant*vector.y, constant*vector.z, constant*vector.w);
}

inline __device__ auto operator*(float2 vector, float constant) -> float2 {
    return make_float2(vector.x*constant, vector.y*constant);
}

inline __device__ auto operator*(float3 vector, float constant) -> float3 {
    return make_float3(vector.x*constant, vector.y*constant, vector.z*constant);
}

inline __device__ auto operator*(float4 vector, float constant) -> float4 {
    return make_float4(vector.x*constant, vector.y*constant, vector.z*constant, vector.w*constant);
}

inline __device__ auto operator*(float constant, float2 vector) -> float2 {
    return make_float2(constant*vector.x, constant*vector.y);
}

inline __device__ auto operator*(float constant, float3 vector) -> float3 {
    return make_float3(constant*vector.x, constant*vector.y, constant*vector.z);
}

inline __device__ auto operator*(float constant, float4 vector) -> float4 {
    return make_float4(constant*vector.x, constant*vector.y, constant*vector.z, constant*vector.w);
}

inline __device__ auto operator*(double2 vector, double constant) -> double2 {
    return make_double2(vector.x*constant, vector.y*constant);
}

inline __device__ auto operator*(double3 vector, double constant) -> double3 {
    return make_double3(vector.x*constant, vector.y*constant, vector.z*constant);
}

inline __device__ auto operator*(double4 vector, double constant) -> double4 {
    return make_double4(vector.x*constant, vector.y*constant, vector.z*constant, vector.w*constant);
}

inline __device__ auto operator*(double constant, double2 vector) -> double2 {
    return make_double2(constant*vector.x, constant*vector.y);
}

inline __device__ auto operator*(double constant, double3 vector) -> double3 {
    return make_double3(constant*vector.x, constant*vector.y, constant*vector.z);
}

inline __device__ auto operator*(double constant, double4 vector) -> double4 {
    return make_double4(constant*vector.x, constant*vector.y, constant*vector.z, constant*vector.w);
}

/** @} */

/**
 * @defgroup vec_scalar_div Scalar-broadcast division (vector / scalar)
 * @brief Divide every component by one scalar. Only the (vector / scalar) order
 *        exists; there is no scalar / vector overload.
 *
 * @param[in] vector   dividend.
 * @param[in] constant scalar divisor.
 * @return the divided vector.
 * @note Precision: the int overloads divide each component (truncating toward
 *       zero). The float/double overloads instead compute scale = 1/constant
 *       once and multiply each component by it; this is NOT bit-identical to
 *       per-component division, because the reciprocal is rounded before the
 *       multiply, adding one rounding step, so a result can differ from true
 *       division by up to one ulp. Callers needing correctly-rounded
 *       per-component division must instead use the vector/vector operator/.
 * @warning A zero @p constant is undefined behavior for the int overloads and
 *          yields inf/nan for float/double; no divisor check is performed.
 * @{
 */
// Divide a vector by a constant.

inline __device__ auto operator/(int2 vector, int constant) -> int2 {
    return make_int2(vector.x/constant, vector.y/constant);
}

inline __device__ auto operator/(int3 vector, int constant) -> int3 {
    return make_int3(vector.x/constant, vector.y/constant, vector.z/constant);
}

inline __device__ auto operator/(int4 vector, int constant) -> int4 {
    return make_int4(vector.x/constant, vector.y/constant, vector.z/constant, vector.w/constant);
}

inline __device__ auto operator/(float2 vector, float constant) -> float2 {
    float scale = 1.0F/constant;
    return vector*scale;
}

inline __device__ auto operator/(float3 vector, float constant) -> float3 {
    float scale = 1.0F/constant;
    return vector*scale;
}

inline __device__ auto operator/(float4 vector, float constant) -> float4 {
    float scale = 1.0F/constant;
    return vector*scale;
}

inline __device__ auto operator/(double2 vector, double constant) -> double2 {
    double scale = 1.0/constant;
    return vector*scale;
}

inline __device__ auto operator/(double3 vector, double constant) -> double3 {
    double scale = 1.0/constant;
    return vector*scale;
}

inline __device__ auto operator/(double4 vector, double constant) -> double4 {
    double scale = 1.0/constant;
    return vector*scale;
}

/** @} */

/**
 * @defgroup vec_scalar_compound Compound scalar scaling (vector *= scalar)
 * @brief In-place scaling of every component by one scalar: vector.c *= constant.
 *
 * @p vector is mutated in place; the operator returns void (no chaining). Only
 * *= is provided for the scalar case — there is no scalar /=, +=, or -=. The
 * vector and scalar element types match per overload; no mixed-precision form.
 *
 * @param[in,out] vector   updated in place with each component scaled.
 * @param[in]     constant scalar multiplier.
 * @{
 */
// *= operator (multiply vector by constant)

inline __device__ void operator*=(int2& vector, int constant) {
    vector.x *= constant; vector.y *= constant;
}

inline __device__ void operator*=(int3& vector, int constant) {
    vector.x *= constant; vector.y *= constant; vector.z *= constant;
}

inline __device__ void operator*=(int4& vector, int constant) {
    vector.x *= constant; vector.y *= constant; vector.z *= constant; vector.w *= constant;
}

inline __device__ void operator*=(float2& vector, float constant) {
    vector.x *= constant; vector.y *= constant;
}

inline __device__ void operator*=(float3& vector, float constant) {
    vector.x *= constant; vector.y *= constant; vector.z *= constant;
}

inline __device__ void operator*=(float4& vector, float constant) {
    vector.x *= constant; vector.y *= constant; vector.z *= constant; vector.w *= constant;
}

inline __device__ void operator*=(double2& vector, double constant) {
    vector.x *= constant; vector.y *= constant;
}

inline __device__ void operator*=(double3& vector, double constant) {
    vector.x *= constant; vector.y *= constant; vector.z *= constant;
}

inline __device__ void operator*=(double4& vector, double constant) {
    vector.x *= constant; vector.y *= constant; vector.z *= constant; vector.w *= constant;
}

/** @} */

/**
 * @brief Euclidean dot (inner) product of two 3-vectors.
 *
 * Returns lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z, the scalar inner product.
 * Accumulated entirely in the element type (float or double) with no
 * intermediate widening, so the float3 overload carries float rounding at every
 * step. Provided only for the 3-component width, only for float/double.
 *
 * @param[in] lhs left operand.
 * @param[in] rhs right operand.
 * @return the scalar inner product lhs·rhs.
 * @note The two products may fuse with the additions under FMA contraction
 *       (nvcc --fmad), shifting the low-order bit relative to strict
 *       multiply-then-add ordering.
 */
// Dot product

inline __device__ auto dot(float3 lhs, float3 rhs) -> float {
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

inline __device__ auto dot(double3 lhs, double3 rhs) -> double {
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

/**
 * @brief Right-handed 3-D vector cross product lhs × rhs.
 *
 * Computes, component-wise:
 *   result.x = lhs.y*rhs.z - lhs.z*rhs.y
 *   result.y = lhs.z*rhs.x - lhs.x*rhs.z
 *   result.z = lhs.x*rhs.y - lhs.y*rhs.x
 * The result is orthogonal to both operands with magnitude |lhs||rhs|·sin θ, and
 * the operation is anti-commutative: cross(a,b) == -cross(b,a).
 *
 * The float4/double4 overloads use only the x/y/z lanes, ignore both inputs' .w,
 * and set result.w = 0; the .w lane is not part of the algebra. They exist so a
 * padded or homogeneous 4-vector can be crossed without unpacking.
 *
 * @param[in] lhs left operand.
 * @param[in] rhs right operand.
 * @return lhs × rhs; result.w is forced to 0 for the 4-component overloads.
 * @note Each component's two products may fuse with its subtraction under FMA
 *       contraction (nvcc --fmad), affecting the last ulp.
 */
// Cross product

inline __device__ auto cross(float3 lhs, float3 rhs) -> float3 {
    return make_float3((lhs.y * rhs.z) - (lhs.z * rhs.y), (lhs.z * rhs.x) - (lhs.x * rhs.z), (lhs.x * rhs.y) - (lhs.y * rhs.x));
}

inline __device__ auto cross(float4 lhs, float4 rhs) -> float4 {
    return make_float4((lhs.y *rhs.z) - (lhs.z * rhs.y), (lhs.z * rhs.x) - (lhs.x * rhs.z), (lhs.x * rhs.y) - (lhs.y * rhs.x), 0.0F);
}

inline __device__ auto cross(double3 lhs, double3 rhs) -> double3 {
    return make_double3((lhs.y * rhs.z) - (lhs.z * rhs.y), (lhs.z * rhs.x) - (lhs.x * rhs.z), (lhs.x * rhs.y) - (lhs.y * rhs.x));
}

inline __device__ auto cross(double4 lhs, double4 rhs) -> double4 {
    return make_double4((lhs.y *rhs.z) - (lhs.z * rhs.y), (lhs.z * rhs.x) - (lhs.x * rhs.z), (lhs.x * rhs.y) - (lhs.y * rhs.x), 0.0);
}

/**
 * @brief Scale a vector to unit length over all of its components.
 *
 * Returns value * rsqrt(Σ_c value.c²), i.e. each component divided by the
 * Euclidean norm taken over every lane of the given width: the 2-wide overloads
 * use x²+y², the 3-wide overloads x²+y²+z², and the 4-wide overloads INCLUDE w
 * (x²+y²+z²+w²). The 4-wide forms are therefore a genuine 4-D normalization, not
 * a 3-vector normalize that leaves the .w lane untouched.
 *
 * @param[in] value vector to normalize.
 * @return @p value scaled to approximately unit length.
 * @note Precision: uses the hardware fast reciprocal square root (rsqrtf for
 *       float, rsqrt for double), so results are approximate (a few ulp) rather
 *       than correctly rounded; the float overloads are the least accurate. The
 *       squared-norm sum is also subject to FMA contraction (nvcc --fmad).
 * @warning No zero-length guard: a zero-length (or denormal-underflowing) input
 *          produces inf/nan components.
 */
// Normalize a vector

inline __device__ auto normalize(float2 value) -> float2 {
    return value*rsqrtf((value.x * value.x) + (value.y * value.y));
}

inline __device__ auto normalize(float3 value) -> float3 {
    return value*rsqrtf((value.x * value.x) + (value.y * value.y) + (value.z * value.z));
}

inline __device__ auto normalize(float4 value) -> float4 {
    return value*rsqrtf((value.x * value.x) + (value.y * value.y) + (value.z * value.z) + (value.w * value.w));
}

inline __device__ auto normalize(double2 value) -> double2 {
    return value*rsqrt((value.x * value.x) + (value.y * value.y));
}

inline __device__ auto normalize(double3 value) -> double3 {
    return value*rsqrt((value.x * value.x) + (value.y * value.y) + (value.z * value.z));
}

inline __device__ auto normalize(double4 value) -> double4 {
    return value*rsqrt((value.x * value.x) + (value.y * value.y) + (value.z * value.z) + (value.w * value.w));
}

/**
 * @brief Drop the fourth component, narrowing a 4-vector to a 3-vector.
 *
 * Returns the x/y/z lanes with value.w discarded — a pure component copy with no
 * arithmetic and no type conversion; the element type is preserved
 * (short4→short3, int4→int3, float4→float3, double4→double3). This is the only
 * function in the file that operates on the short vector types.
 *
 * @param[in] value 4-component source.
 * @return the x/y/z lanes as a 3-component vector.
 */
// Strip off the fourth component of a vector.

inline __device__ auto trimTo3(short4 value) -> short3 {
    return make_short3(value.x, value.y, value.z);
}

inline __device__ auto trimTo3(int4 value) -> int3 {
    return make_int3(value.x, value.y, value.z);
}

inline __device__ auto trimTo3(float4 value) -> float3 {
    return make_float3(value.x, value.y, value.z);
}

inline __device__ auto trimTo3(double4 value) -> double3 {
    return make_double3(value.x, value.y, value.z);
}
