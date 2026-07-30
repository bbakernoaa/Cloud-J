# Tasks: Optional PCR Tridiagonal Solver Integration

**Input**: Design documents from `specs/005-pcr-optional-solver/`

**Prerequisites**: plan.md (required), spec.md (required for user stories)

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

---

## Phase 1: Setup

**Purpose**: Conditionally declare the PCR compiler flag.

- [x] T001 Add `CLOUDJ_USE_PCR` option to `CMakeLists.txt` and conditionally add compilation definition `CLOUDJ_USE_PCR` when enabled.

---

## Phase 2: Foundational

**Purpose**: Define the block-tridiagonal PCR solver mathematics.

- [x] T002 Implement the block-tridiagonal Parallel Cyclic Reduction (PCR) algorithm `solve_pcr` in `include/cloudj/radiative_solver.hpp` wrapped inside `#if defined(CLOUDJ_USE_PCR)` blocks.

---

## Phase 3: User Story 2 - GPU Vertical Parallelization (Priority: P1) 🎯 MVP

**Goal**: Redirect block-tridiagonal execution to PCR if enabled.

- [x] T003 Conditionally redirect standard `BLKSLV` execution inside `include/cloudj/radiative_solver.hpp` to the PCR solver backend if `CLOUDJ_USE_PCR` compiler flag is active.

---

## Phase 4: Polish & Cross-Cutting Concerns

**Purpose**: Validation and compiling checks.

- [x] T004 Compile standard build (`-DCLOUDJ_USE_PCR=OFF`) and verify default Thomas-Feautrier compilation remains clean.
- [x] T005 Compile PCR-enabled build (`-DCLOUDJ_USE_PCR=ON`) and verify build succeeds.
- [x] T006 Run the 500,000 fuzzer stress checks to guarantee exact double-precision numerical parity remains intact across both execution paths.

---

## Dependencies & Execution Order

1. **Phase 1 (Setup)**: Prerequisite for all tasks.
2. **Phase 2 (Foundational)**: Core PCR solver implementation.
3. **Phase 3 (User Story 2)**: Conditional execution redirect.
4. **Phase 4 (Polish)**: Validation.
