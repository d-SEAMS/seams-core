<!-- https://docs.dseams.info/api/class_seams_Error.md -->

# class seams::Error {#class-seams-error}

## Overview {#overview}

Failure from the C API, or from a rejected argument. [More…](#details-classseams-1-1-error)

## Detailed Documentation {#detailed-documentation}

Failure from the C API, or from a rejected argument.

`what()` is the text from `seams_handle_error`. `code()` is the status integer. The object does not own a handle.

### Construction {#construction}

Build an error from a status and a message.

**Parameters:**

| status | Status code. Zero is not thrown by check. |
| --- | --- |
| message | Detail text. An empty string is allowed. |

### Methods {#methods}

[Status](generated_targets.md#doxid-group-seams-cpp-1ga66e498496a8ac061410534a9f0171779) as the `Status` enumeration when the code is known.

**Returns:**

The matching enumerator.

[Status](generated_targets.md#doxid-group-seams-cpp-1ga66e498496a8ac061410534a9f0171779) integer.

Same values as `seams_status`.

**Returns:**

The code stored at construction.
