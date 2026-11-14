# Rust Module Style Guide
This doc provides guidelines to follow when modifying or adding to this codebase. These are meant to be guidelines, not necessarily hard and fast rules. As always, use your best judgement. The goal of style, structure, comments, and documentation is to eliminate confusion.

## Module Structure
Each rust module should follow this general structure

```rust
//! Description of the module

// ------------------
// External Libraries
// ------------------

// ------------------
// Internal Libraries
// ------------------

// -------
// Structs
// -------

// -----
// Enums
// -----

// ------
// Traits
// ------

// ---------
// Constants
// ---------

// ---------
// Functions
// ---------

// ----------
// Unit Tests
// ----------

```

Structs should follow this structure

```rust
/// Brief struct title
///
/// Detailed struct description
///
/// # Examples
/// ```rust
/// ```
///
/// # References
/// - [Source](link)
#[derive(attributes)]
pub struct ExampleStruct {
    /// Variable description
    pub example_var_1: type,

    /// Variable description
    pub example_var_2: type,
}
```

Enums should follow this structure

```rust
/// Brief enum title
///
/// Detailed enum description
///
/// # Examples
/// ```rust
/// ```
///
/// # References
/// - [Source](link)
#[derive(attributes)]
pub enum ExampleEnum {
    /// Enum description
    #[default]
    Enum1,

    /// Enum description
    Enum2,
}
```

Public enums are marked `#[non_exhaustive]` after the derive, so adding a variant is not a breaking change

Constants should follow this structure

```rust
/// Brief constant title
///
/// Detailed constant description
///
/// # Examples
/// ```rust
/// ```
///
/// # References
/// - [Source](link)
const CONSTANT1: type = value;
```

Functions should follow this structure

```rust
/// Brief function title
///
/// Detailed function description
///
/// # Arguments
/// * `arg1` - Argument description
///
/// # Returns
/// * Return description
///
/// # Panics
/// * Panic descriptions
///
/// # Errors
/// * Error descriptions
///
/// # Safety
/// * Safety descriptions
///
/// # Examples
/// ```rust
/// ```
///
/// # References
/// - [Source](link)
pub fn example_function(arg1: type) -> (return1) {
    // Well commented
}
```

Note that for all non-`pub`, the `# Examples` section is unnecessary.

## Special Characters
Avoid including special characters in code and docstrings.