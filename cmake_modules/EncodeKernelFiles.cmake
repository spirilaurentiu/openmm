file(GLOB KERNEL_FILES ${KERNEL_SOURCE_DIR}/kernels/*.${KERNEL_FILE_EXTENSION})
set(KERNEL_FILE_DECLARATIONS)
configure_file(${KERNEL_SOURCE_DIR}/${KERNEL_SOURCE_CLASS}.cpp.in
               ${KERNELS_CPP})
# Determine file extension length
string(LENGTH ${KERNEL_FILE_EXTENSION} extension_length)
# add one space for the dot
math(EXPR extension_length ${extension_length}+1)
# Locate a Python interpreter used to strip comments from each kernel before
# encoding. Kernel comments would otherwise be baked verbatim into the generated
# string literal, ship in the .so, and be re-lexed by NVRTC on every cold JIT
# (they also enter the SHA1 the JIT cache is keyed on). Stripping preserves line
# numbers (see strip_comments.py). If no interpreter is found we fall back to
# encoding the raw file, so the build never breaks over this.
find_program(ENCODE_KERNELS_PYTHON NAMES python3 python)
set(ENCODE_KERNELS_STRIP_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/strip_comments.py)
foreach(file ${KERNEL_FILES})
  # Load the file contents (comments stripped when a Python interpreter is
  # available, raw otherwise) and process it.
  set(kernel_stripped FALSE)
  if(ENCODE_KERNELS_PYTHON AND EXISTS ${ENCODE_KERNELS_STRIP_SCRIPT})
    execute_process(
      COMMAND ${ENCODE_KERNELS_PYTHON} ${ENCODE_KERNELS_STRIP_SCRIPT} ${file}
      OUTPUT_VARIABLE file_content
      RESULT_VARIABLE kernel_strip_result)
    if(kernel_strip_result EQUAL 0)
      set(kernel_stripped TRUE)
    endif()
  endif()
  if(NOT kernel_stripped)
    file(STRINGS ${file} file_content NEWLINE_CONSUME)
  endif()
  # Replace all backslashes by double backslashes as they are being put in a C
  # string. Be careful not to replace the backslash before a semicolon as that
  # is the CMAKE internal escaping of a semicolon to prevent it from acting as a
  # list separator.
  string(REGEX REPLACE "\\\\([^;])" "\\\\\\\\\\1" file_content
                       "${file_content}")
  # Escape double quotes as being put in a C string.
  string(REPLACE "\"" "\\\"" file_content "${file_content}")
  # Split in separate C strings for each line.
  string(REPLACE "\n" "\\n\"\n\"" file_content "${file_content}")

  # Determine a name for the variable that will contain this file's contents
  file(RELATIVE_PATH filename ${KERNEL_SOURCE_DIR}/kernels ${file})
  string(LENGTH ${filename} filename_length)
  math(EXPR filename_length ${filename_length}-${extension_length})
  string(SUBSTRING ${filename} 0 ${filename_length} variable_name)

  # Record the variable declaration and definition.
  set(KERNEL_FILE_DECLARATIONS
      ${KERNEL_FILE_DECLARATIONS}static\ const\ std::string\ ${variable_name};\n
  )
  file(
    APPEND ${KERNELS_CPP}
    const\ string\ ${KERNEL_SOURCE_CLASS}::${variable_name}\ =\ \"${file_content}\"\;\n
  )
endforeach(file)
configure_file(${KERNEL_SOURCE_DIR}/${KERNEL_SOURCE_CLASS}.h.in ${KERNELS_H})
