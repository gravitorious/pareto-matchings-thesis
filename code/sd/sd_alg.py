#serial dictatorship
#complexity: O(n!)
import csv
import itertools
from pathlib import Path

def read_preferences(filename):
    preferences = {}
    with open(filename, newline="") as prefcsv:
        reader = csv.DictReader(prefcsv)
        for row in reader:
            applicant = int(row["applicant"])
            prefs = []
            for key, value in row.items():
                if key.startswith("pref") and value != "":
                    prefs.append(int(value))
            preferences[applicant] = prefs
    return preferences

def matching_to_tuple(matching, applicants):
    return{
        a: matching.get(a, -1)
        for a in sorted(applicants)
    }

def serial_dictatorship(preferences):
    matchings = {}
    m = len(preferences) #applicants (number of rows)
    applicants = list(range(1, m + 1))
    for perm in itertools.permutations(applicants):
        assigned_houses = set()
        matching = {}
        for applicant in perm:
            for house in preferences[applicant]:
                if house not in assigned_houses:
                    matching[applicant] = house
                    assigned_houses.add(house)
                    break
        matching = matching_to_tuple(matching, applicants)
        matchings[perm] = matching
    return matchings

def all_solutions(matchings, instance):
    #matchings = {perm, dict with matchings}
    n = len(next(iter(matchings.values())))  #applicants
    path = Path("results/" + instance + "/" + instance + "_sd_pareto_matchings.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as sdcsv:
        writer = csv.writer(sdcsv)
        writer.writerow(["order", "matching"])
        for perm, matching in matchings.items():
            order_str = ",".join(map(str, perm))
            alloc_list = [matching.get(a, -1) for a in perm]
            alloc_str = ",".join(map(str, alloc_list))
            writer.writerow([order_str, alloc_str])

def get_unique_matchings(matchings, m):
    unique = set()
    for matching in matchings.values():
        alloc_tuple = tuple(
            matching.get(a, -1) for a in range(1, m + 1)
        )
        unique.add(alloc_tuple)
    return sorted(unique)

def write_unique_matchings(matchings, m, instance):
    path = Path("results/" + instance + "/" + instance + "_unique_pareto_matchings.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    unique = get_unique_matchings(matchings, m)
    with open(path, "w", newline="") as uniquecsv:
        writer = csv.writer(uniquecsv)
        writer.writerow([f"A{i}" for i in range(1, m + 1)])
        for alloc in sorted(unique):
            writer.writerow(alloc)
    return unique

def write_pom_vectors(unique_matchings, preferences, instance):
    col = []
    #create Xij columns
    for applicant in sorted(preferences.keys()):
        for house in sorted(set(preferences[applicant])):
            col.append((applicant, house))
    path = Path("results/" + instance + "/" + instance + "_pom_vectors.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as vectorcsv:
        writer = csv.writer(vectorcsv)
        #header
        header = [f"X{applicant}{house}" for applicant, house in col]
        writer.writerow(header)
        #rows
        for matching in unique_matchings:
            row = []
            for applicant, house in col:
                assigned_house = matching[applicant - 1]
                if assigned_house == house:
                    row.append(1)
                else:
                    row.append("")
            writer.writerow(row)

def compute_matching_statistics(unique_matchings, preferences):
    stats = {}
    #number of unique matchings
    stats["num_matchings"] = len(unique_matchings)
    #size of each matching = number of assigned applicants
    matching_sizes = [
        sum(1 for house in matching if house != -1)
        for matching in unique_matchings
    ]
    min_matching = min(matching_sizes)
    max_matching = max(matching_sizes)
    stats["min_matching"] = min_matching
    stats["max_matching"] = max_matching
    #matching numbers with minimum / maximum size
    stats["min_size_matchings"] = [
        i + 1
        for i, size in enumerate(matching_sizes)
        if size == min_matching
    ]
    stats["max_size_matchings"] = [
        i + 1
        for i, size in enumerate(matching_sizes)
        if size == max_matching
    ]
    #variables Xij that are never equal to 1 in any unique matching
    not_used = []
    for applicant in sorted(preferences.keys()):
        for house in sorted(set(preferences[applicant])):
            used = False
            for matching in unique_matchings:
                if matching[applicant - 1] == house:
                    used = True
                    break
            if not used:
                not_used.append((applicant, house))
    stats["not_used"] = not_used
    #variables Xij that are equal to 1 in every unique matching
    always_used = []
    for applicant in sorted(preferences.keys()):
        for house in sorted(set(preferences[applicant])):
            used_in_all = True
            for matching in unique_matchings:
                if matching[applicant - 1] != house:
                    used_in_all = False
                    break
            if used_in_all:
                always_used.append((applicant, house))
    stats["always_used"] = always_used
    return stats

def write_matching_statistics_txt(stats, instance):
    path = Path("results/" + instance + "/" + instance + "_matching_statistics.txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as stattxt:
        stattxt.write("STATS\n")
        stattxt.write(f"Number of unique matchings: {stats['num_matchings']}\n")
        stattxt.write(
            f"Minimum matching size: {stats['min_matching']} "
            f"(matchings: {', '.join(map(str, stats['min_size_matchings']))})\n"
        )
        stattxt.write(
            f"Maximum matching size: {stats['max_matching']} "
            f"(matchings: {', '.join(map(str, stats['max_size_matchings']))})\n"
        )
        for applicant, house in stats["not_used"]:
            stattxt.write(f"X{applicant}{house} is never equal to 1 in any matching.\n")
        for applicant, house in stats["always_used"]:
            stattxt.write(f"X{applicant}{house} is equal to 1 in every matching.\n")

def sd_alg(instance):
    input = "data/" + instance + ".csv"
    preferences = read_preferences(input)
    matchings = serial_dictatorship(preferences)
    m = len(preferences) #applicants
    all_solutions(matchings, instance)
    write_unique_matchings(matchings, len(preferences), instance)
    unique_matchings = get_unique_matchings(matchings, m)
    write_pom_vectors(unique_matchings, preferences, instance)
    stats = compute_matching_statistics(unique_matchings, preferences)
    write_matching_statistics_txt(stats, instance)
