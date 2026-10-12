#include "common.h"
#include "grouping.h"

#include <algorithm>
#include <numeric>

using namespace Rcpp;

namespace tsir {

std::vector<int>
group_by_linesweep(
    const std::vector<int64_t> &start, const std::vector<int64_t> &end, const std::vector<double> &time)
{
    const size_t n_in = start.size();
    std::vector<int> group_of(n_in, -1);
    if (n_in == 0) {
        return group_of;
    }
    /* Merge overlapping same-time haplotypes (sorted by time, then start). */
    std::vector<size_t> order(n_in);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (time[a] != time[b]) {
            return time[a] < time[b];
        }
        return start[a] < start[b];
    });
    std::vector<int64_t> new_start, new_end;
    std::vector<double> new_time;
    std::vector<std::vector<size_t>> members;
    size_t i = 0;
    while (i < n_in) {
        size_t j = i + 1;
        std::vector<size_t> group{ order[i] };
        int64_t max_right = end[order[i]];
        while (j < n_in && time[order[j]] == time[order[i]] && start[order[j]] < max_right) {
            max_right = std::max(max_right, end[order[j]]);
            group.push_back(order[j]);
            j++;
        }
        new_start.push_back(start[order[i]]);
        new_end.push_back(max_right);
        new_time.push_back(time[order[i]]);
        members.push_back(std::move(group));
        i = j;
    }
    const size_t n = new_time.size();

    /* Linesweep over start/end events, ends before starts at equal positions. */
    std::vector<size_t> events(2 * n);
    std::iota(events.begin(), events.end(), 0);
    auto event_pos = [&](size_t e) { return e < n ? new_start[e] : new_end[e - n]; };
    auto event_type = [&](size_t e) { return e < n ? 1 : 0; };
    std::stable_sort(events.begin(), events.end(), [&](size_t a, size_t b) {
        if (event_pos(a) != event_pos(b)) {
            return event_pos(a) < event_pos(b);
        }
        return event_type(a) < event_type(b);
    });
    std::vector<std::vector<int>> children(n);
    std::vector<int> incoming(n, 0);
    std::vector<int> active;
    for (size_t e : events) {
        int index = static_cast<int>(e < n ? e : e - n);
        if (event_type(e) == 1) {
            double e_time = new_time[index];
            for (int a : active) {
                if (new_time[a] > e_time) {
                    incoming[index]++;
                    children[a].push_back(index);
                } else if (new_time[a] < e_time) {
                    incoming[a]++;
                    children[index].push_back(a);
                }
            }
            active.push_back(index);
        } else {
            auto it = std::find(active.begin(), active.end(), index);
            if (it != active.end()) {
                active.erase(it);
            }
        }
    }

    /* Topological sort, one level per group. */
    std::vector<int> group_id(n, -1);
    int current = 0;
    std::vector<int> ready;
    while (true) {
        ready.clear();
        for (size_t u = 0; u < n; u++) {
            if (incoming[u] == 0) {
                ready.push_back(static_cast<int>(u));
            }
        }
        if (ready.empty()) {
            break;
        }
        for (int u : ready) {
            incoming[u] = -1;
            for (int c : children[u]) {
                incoming[c]--;
            }
            group_id[u] = current;
        }
        current++;
    }
    for (size_t u = 0; u < n; u++) {
        if (group_id[u] == -1) {
            throw std::runtime_error(
                "Erroneous cycle in ancestor dependancies, this is often caused by too "
                "many unique site times. This fixed by discretising the site times, for "
                "example rounding times to the nearest 0.1.");
        }
        for (size_t m : members[u]) {
            group_of[m] = group_id[u];
        }
    }
    return group_of;
}

} // namespace tsir

// [[Rcpp::export(.group_by_linesweep)]]
IntegerVector
group_by_linesweep_r(NumericVector start, NumericVector end, NumericVector time)
{
    if (start.size() != end.size() || start.size() != time.size()) {
        stop("start, end and time must have the same length");
    }
    std::vector<int64_t> s(start.begin(), start.end()), e(end.begin(), end.end());
    std::vector<double> t(time.begin(), time.end());
    try {
        std::vector<int> g = tsir::group_by_linesweep(s, e, t);
        return IntegerVector(g.begin(), g.end());
    } catch (std::exception &ex) {
        stop(ex.what());
    }
}
