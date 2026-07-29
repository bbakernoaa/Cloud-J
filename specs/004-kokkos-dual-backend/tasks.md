# Tasks: Kokkos Dual-Backend Integration

**Input**: Design documents from `specs/004-kokkos-dual-backend/`

**Prerequisites**: plan.md (required), spec.md (required for user stories)

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

---

## Phase 1: Setup

**Purpose**: Conditionally locate Kokkos during CMake configuring.

- [x] T001 Add `CLOUDJ_USE_KOKKOS` option to `CMakeLists.txt` and conditionally load Kokkos and KokkosKernels packages if enabled.

---

## Phase 2: Foundational

**Purpose**: Define the conditional Kokkos View data structures.

- [x] T002 Create `include/cloudj/kokkos_backend.hpp` and declare `Kokkos::View` layouts mapping to standard Cloud-J profiles and rates matrices inside `#if defined(CLOUDJ_USE_KOKKOS)` blocks.

---

## Phase 3: User Story 2 - Exascale Performance Portability (Priority: P1) 🎯 MVP

**Goal**: Implement Kokkos parallel dispatch for performance-portable solves.

- [x] T003 Implement portable execution wrapper in `include/cloudj/kokkos_backend.hpp` to offload column calculations using `Kokkos::parallel_for`.
- [x] T004 Wrap standard photolysis loops inside `calculate_photolysis_rates` in `include/cloudj/cloudj.hpp` with this execution policy wrapper.

---

## Phase 4: User Story 3 - National Lab Solver Acceleration via KokkosKernels (Priority: P2)

**Goal**: Map tridiagonal block steps to KokkosKernels.

- [x] T005 Map matrix solver steps inside `kokkos_backend.hpp` to `KokkosBatched::LU` and `KokkosBatched::Trsv` when KokkosKernels is linked.

---

## Phase 5: Polish & Cross-Cutting Concerns

**Purpose**: Validation and compiling checks.

- [x] T006 Compile standard build (`-DCLOUDJ_USE_KOKKOS=OFF`) and verify zero dependency compilation succeeds.
- [x] T007 Compile Kokkos-enabled build (`-DCLOUDJ_USE_KOKKOS=ON`) and verify linking target consistency.
- [x] T008 Run the 500,000 fuzzer stress checks to guarantee exact double-precision numerical parity remains intact across both execution paths.

---

## Dependencies & Execution Order

1. **Phase 1 (Setup)**: Prerequisite for all tasks.
2. **Phase 2 (Foundational)**: View structures are needed by the execution policies.
3. **Phase 3 (User Story 2)**: Core Kokkos offloading loops.
4. **Phase 4 (User Story 3)**: Acceleration via KokkosKernels.
5. **Phase 5 (Polish)**: Validation.
