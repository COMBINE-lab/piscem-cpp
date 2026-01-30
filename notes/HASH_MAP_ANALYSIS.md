# Hash Map Replacement Analysis: phmap::flat_hash_map → ankerl::unordered_dense::map

## Executive Summary

This document analyzes the feasibility and implications of replacing `phmap::flat_hash_map` with `ankerl::unordered_dense::map` in the piscem-cpp codebase, with a focus on the critical `hit_map` in the mapping hot path.

## Current State

### Library Availability
- ✅ **ankerl::unordered_dense v4.5.0** is already included in `include/unordered_dense.h`
- ✅ **parallel_hashmap (phmap)** is currently used throughout the codebase
- 📝 One commented-out usage of `ankerl::unordered_dense::map` exists in `streaming_query.hpp` (line 36)

### Primary Hash Map: hit_map

**Location**: `include/mapping/utils.hpp` line 947 in `mapping_cache_info` template class

**Type**: `phmap::flat_hash_map<uint32_t, sketch_hit_info_t>`

**Purpose**: Core mapping structure that stores sketch hit information indexed by transcript/target ID

### All phmap::flat_hash_map Usage

| Location | Type | Purpose | Hot Path? |
|----------|------|---------|-----------|
| `mapping/utils.hpp` | `flat_hash_map<uint32_t, sketch_hit_info_t>` | **hit_map** - primary mapping cache | ✅ YES |
| `mapping/utils.hpp` | `flat_hash_map<uint64_t, uint32_t>` | unmapped_bc_map - barcode tracking | ⚠️ Moderate |
| `poison_table.hpp` | `flat_hash_map<uint64_t, uint64_t, RobinHoodHash>` | poison_map_t - poison k-mer storage | ❌ No (build-time) |
| `index_evaluator.cpp` | `flat_hash_map<uint64_t, uint64_t>` | freq_map - frequency counting | ❌ No (analysis tool) |
| `build_contig_table.cpp` | `flat_hash_map<uint64_t, rank_count>` | id_to_rank - contig table building | ❌ No (build-time) |
| `build_contig_table.cpp` | `flat_hash_map<vector<...>, rank_offset>` | ec_id_map - equivalence class mapping | ❌ No (build-time) |

## API Compatibility Analysis

### Operations Used on hit_map

| Operation | phmap::flat_hash_map | ankerl::unordered_dense::map | Compatible? |
|-----------|---------------------|------------------------------|-------------|
| `reserve(size)` | ✅ Yes | ✅ Yes | ✅ YES |
| `clear()` | ✅ Yes | ✅ Yes | ✅ YES |
| `operator[key]` | ✅ Yes (insert/access) | ✅ Yes (insert/access) | ✅ YES |
| `find(key)` | ✅ Yes | ✅ Yes | ✅ YES |
| `end()` | ✅ Yes | ✅ Yes | ✅ YES |
| `empty()` | ✅ Yes | ✅ Yes | ✅ YES |
| `size()` | ✅ Yes | ✅ Yes | ✅ YES |
| Range-based for | ✅ Yes | ✅ Yes | ✅ YES |
| Iterator deref (`kv.first`, `kv.second`) | ✅ Yes | ✅ Yes | ✅ YES |

**Conclusion**: ✅ **100% API compatible** - no code changes needed beyond the type declaration

### Custom Hash Function Support

The `poison_map_t` uses a custom hash function:
```cpp
phmap::flat_hash_map<uint64_t, uint64_t, sshash::RobinHoodHash>
```

**ankerl::unordered_dense** also supports custom hash functions:
```cpp
ankerl::unordered_dense::map<uint64_t, uint64_t, sshash::RobinHoodHash>
```

✅ **Compatible** - custom hash template parameter is supported

## Performance Characteristics

### phmap::flat_hash_map

**Strengths:**
- Open addressing with flat storage (good cache locality)
- Linear probing variant
- Low memory overhead
- SIMD optimizations for probe sequence
- Mature, battle-tested library

**Weaknesses:**
- Can suffer from clustering under certain hash distributions
- Delete operations can degrade performance over time
- Not as cache-friendly as some newer designs

### ankerl::unordered_dense::map

**Strengths:**
- **Robin Hood backward shift deletion** - maintains good performance after deletes
- **Extremely cache-friendly**: separate index + data storage
  - Small index array (1 byte per bucket)
  - Dense data vector (contiguous storage)
- **Fast iteration** - data is contiguous, not scattered
- **Lower memory overhead** - typically 30-50% less than std::unordered_map
- **Better worst-case**: bounded probe sequence length
- **Modern C++17 design** with excellent performance characteristics

**Weaknesses:**
- Slightly more complex deletion (backward shift)
- Less battle-tested than phmap in bioinformatics workloads
- May have different performance under extreme load factors

### Expected Performance Impact

Based on the design characteristics and published benchmarks:

#### For hit_map (the critical path):

**Likely Improvements** (⬆️):
- ⬆️ **Iteration speed**: 20-50% faster due to contiguous data storage
  - Used in: mapping result processing (lines 1556+)
  - Impact: Moderate-High
  
- ⬆️ **Memory usage**: 20-40% reduction
  - Less memory → better cache utilization → faster overall
  - Impact: Moderate
  
- ⬆️ **Clear + rebuild patterns**: 10-30% faster
  - hit_map is cleared and rebuilt for every read
  - Contiguous storage means less fragmentation
  - Impact: High (this is the hot path!)

**Potential Regressions** (⬇️):
- ⬇️ **Single insert/lookup**: 0-10% slower in some cases
  - Two-level lookup (index → data) vs. direct
  - Impact: Low (typically offset by iteration gains)

**Overall Expected Impact**: 
- 🎯 **3-12% speedup** in mapping throughput
- 🎯 **15-30% reduction** in memory footprint for mapping cache
- 🎯 **Better scalability** with varying hit counts

## Implementation Requirements

### Minimal Changes Required

#### 1. Update include directive (6 files)
```cpp
// Change from:
#include "parallel_hashmap/phmap.h"

// Change to:
#include "unordered_dense.h"
```

#### 2. Update type declarations (8-10 locations)
```cpp
// Change from:
phmap::flat_hash_map<K, V>

// Change to:
ankerl::unordered_dense::map<K, V>
```

#### 3. Update namespace aliases (optional, for brevity)
```cpp
namespace dense = ankerl::unordered_dense;
// Then use: dense::map<K, V>
```

### Files Requiring Changes

**High Priority (hot path):**
1. ✅ `include/mapping/utils.hpp` - hit_map, unmapped_bc_map

**Medium Priority:**
2. ⚠️ `include/poison_table.hpp` - poison_map_t (with custom hash)

**Low Priority (not hot path, but could standardize):**
3. `src/index_evaluator.cpp` - freq_map
4. `src/build_contig_table.cpp` - id_to_rank, ec_id_map
5. `src/pesc_sc.cpp` - if using phmap

### Compatibility Notes

- ✅ No API changes needed - drop-in replacement
- ✅ No algorithm changes needed
- ✅ Custom hash functions work identically
- ✅ All operations used in codebase are supported
- ⚠️ Need to test with RobinHoodHash custom hash function
- ⚠️ May need load factor tuning for optimal performance

## Other Optimization Opportunities

### 1. Standard Library Maps → ankerl::unordered_dense

**Currently using std::unordered_map:**
- `src/pesc_bulk.cpp` - parameter maps (string → string)
- `src/pesc_sc.cpp` - freq_map (uint32_t → size_t)

**Benefit**: 20-40% speedup, 30-50% memory reduction

**Priority**: Low (not in hot path, but easy win)

### 2. Concurrent Access Patterns

**Current**: `boost::concurrent_flat_map` in streaming_query.hpp

**Already noted**: Comment suggests ankerl::unordered_dense was considered

**Consideration**: ankerl::unordered_dense is NOT thread-safe by default
- If thread-safe access needed, keep boost::concurrent_flat_map
- OR use external synchronization with ankerl::unordered_dense

### 3. Small Vector Optimization (already done!)

✅ **Already optimized**: `itlib::small_vector<uint32_t, 255>` for ambiguous_hit_indices
- This is excellent - avoids heap allocation for small cases

### 4. Hash Set Usage

**Current**: `phmap::flat_hash_set<uint64_t>` for observed_ecs

**Could migrate** to: `ankerl::unordered_dense::set<uint64_t>`

**Benefit**: Similar to map - better iteration, lower memory

## Testing Strategy

### 1. Correctness Testing
- [ ] Run existing test suite with replacement
- [ ] Verify mapping results are identical
- [ ] Test with various input sizes
- [ ] Test edge cases (empty maps, single entries, max capacity)

### 2. Performance Testing
- [ ] Benchmark mapping throughput on representative datasets
- [ ] Measure memory usage (RSS, peak allocation)
- [ ] Profile hot path with perf/vtune
- [ ] Compare iteration vs. lookup performance

### 3. Regression Testing
- [ ] Test with custom hash functions (poison_table)
- [ ] Verify reserve() behavior is optimal
- [ ] Test clear() + rebuild pattern performance

## Recommended Implementation Approach

### Phase 1: Proof of Concept (Low Risk)
1. Create feature branch
2. Replace hit_map in mapping/utils.hpp
3. Run tests and benchmarks
4. Gather performance data

### Phase 2: Extended Implementation (Medium Risk)
1. Replace unmapped_bc_map
2. Replace other mapping-related maps
3. Re-run benchmarks

### Phase 3: Full Migration (Optional)
1. Replace remaining phmap usage
2. Consider removing phmap dependency
3. Standardize on ankerl::unordered_dense

## Risks and Mitigations

### Risk 1: Performance Regression
**Mitigation**: Benchmark before merge; keep easy rollback

### Risk 2: Unexpected API Differences
**Mitigation**: Comprehensive testing; API review shows compatibility

### Risk 3: Custom Hash Compatibility
**Mitigation**: Test poison_table separately; validate RobinHoodHash works

### Risk 4: Memory Usage Patterns
**Mitigation**: Profile actual memory usage; may need load factor tuning

## Conclusion

### Feasibility: ✅ HIGH
- Library already available
- API fully compatible
- Minimal code changes required

### Expected Benefit: ✅ MODERATE TO HIGH
- 3-12% throughput improvement likely
- 15-30% memory reduction expected
- Better performance under varying workloads

### Risk Level: ✅ LOW
- Easy rollback (single type change)
- No algorithm changes
- Extensive test coverage available

### Recommendation: ✅ **PROCEED WITH PHASE 1**

Implement hit_map replacement as proof of concept, benchmark thoroughly, and make data-driven decision for full migration.

## References

- ankerl::unordered_dense: https://github.com/martinus/unordered_dense
- Benchmarks: https://martin.ankerl.com/2022/08/27/hashmap-bench-01/
- Design rationale: Robin Hood hashing with backward shift deletion
- Version in repo: v4.5.0 (latest stable)
