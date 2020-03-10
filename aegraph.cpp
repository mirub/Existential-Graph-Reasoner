// Copyright 2019 Banu Miruna-Elena, Descultu Cristian Petrisor
#include <vector>
#include <stack>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }


    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

void AEGraph::possible_double_cuts_aux(
    std::vector<std::vector<int>> &paths, std::vector<int>
        &current) const {
    // auxiliary function that modifies the paths & current
    // vectors
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    if (node->subgraphs.size()) {
        for (int i = 0; (unsigned) i < node->subgraphs.size(); ++i) {
            // adds the cut to the current positions vector
            current.push_back({i});
            if (node->subgraphs[i].subgraphs.size() == 1 &&
                !node->subgraphs[i].atoms.size() && current.size()) {
            // is a double cut positions and adds the current vector
            // to the paths vector
                paths.push_back(current);
            }
            node->subgraphs[i].possible_double_cuts_aux(paths, current);
        }
        // deletes the last cut
        current.pop_back();
    } else {
        if (current.size()) {
            // deletes the last cut
            current.pop_back();
        }
    }
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // 10p
    // creates the vectors needed for the auxiliary function and
    // modifies them upon calling
    std::vector<std::vector<int>> paths;
    std::vector<int> current;
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    node->possible_double_cuts_aux(paths, current);
    // returns all the positions
    return paths;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // 10p
    // creates a pointer to the current graph
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    for (int i = 0; (unsigned) i < where.size() - 1; ++i) {
        // places the pointer to the position above the
        // double cut
        node = &node->subgraphs[where[i]];
    }
    // auxiliary pointer towards the node in the cuts
    AEGraph *aux = &node->subgraphs[where[where.size() - 1]].subgraphs[0];
    for (int i = 0; (unsigned) i < aux->subgraphs.size(); ++i) {
        // adds the subgraphs of the 'aux' node to the node
        // 2 cuts above
        node->subgraphs.push_back(aux->subgraphs[i]);
    }
    for (int i = 0; (unsigned) i < aux->atoms.size(); ++i) {
        // adds the atoms of the 'aux' node to the node
        // 2 cuts above
        node->atoms.push_back(aux->atoms[i]);
    }
    // deletes the double cut subgraph
    node->subgraphs.erase(node->subgraphs.begin() + where[where.size() - 1]);
    return node_value;
}

void AEGraph::possible_erasures_aux(
    std::vector<std::vector<int>> &paths, std::vector<int>
        &current, int level) const {
    // auxiliary function that modifies the paths & current
    // vectors
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    for (int i = 0; (unsigned) i < (node->subgraphs.size() +
        node->atoms.size()); ++i) {
        // adds the cut to the current vector
        current.push_back({i});
        // if the current level is an even one
        if (level % 2 != 0) {
            // adds the current vector to the paths vector
            paths.push_back(current);
            if (level != -1 && (!node->subgraphs.size() &&
                node->atoms.size() == 1) && paths.size()) {
                // deletes the last entry
                paths.pop_back();
            }
            if (level != -1 && (node->subgraphs.size() == 1
                && !node->atoms.size()) && paths.size()) {
                // deletes the last entry
                paths.pop_back();
            }
        }
        if ((unsigned) i < subgraphs.size()) {
            // calls recursivity for the subgraphs
            node->subgraphs[i].possible_erasures_aux(paths, current, level + 1);
        } else {
            // deletes the last cut
            current.pop_back();
        }
    }
    current.pop_back();
}

std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    // 10p
    // creates the vectors needed for the auxiliary function and
    // modifies them upon calling
    std::vector<std::vector<int>> paths;
    std::vector<int> current;
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    // returns all the positions
    node->possible_erasures_aux(paths, current, level);
    return paths;
}

AEGraph AEGraph::erase(std::vector<int> where) const {
    // 10p
    // creates a pointer to the current graph
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    for (int i = 0; (unsigned) i < where.size() - 1; ++i) {
        // places the pointer above the erase position
        node = &node->subgraphs[where[i]];
    }
    if (node->subgraphs.size() <= (unsigned) where[where.size() - 1]) {
        // if the following cut is on a subgraph, erase it
        node->atoms.erase(node->atoms.begin() +
            (where[where.size() - 1] - node->subgraphs.size()));
    } else{
        // if the following cut is on an atom, erase it
        node->subgraphs.erase(node->subgraphs.begin() +
            where[where.size() - 1]);
    }
    return node_value;
}

void AEGraph::possible_deiterations_aux(
        std::vector<std::vector<int>> &paths) const {
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    std::vector<std::vector<int>> way;
    for (int i = 0; (unsigned) i < node->subgraphs.size(); ++i) {
        for (int j = 0; (unsigned) j < node->subgraphs.size(); ++j) {
            // searched for the current graph in the others
            if (i != j &&
                node->subgraphs[j].contains(node->subgraphs[i])) {
                // if found, gets the ways to it
                way = node->subgraphs[j].get_paths_to(node->subgraphs[i]);
                for (auto& v : way) {
                    // adds the previous cut to every way returned by the
                    // get_paths_to function
                    v.insert(v.begin(), j);
                    // adds the way vector to the paths vector
                    paths.push_back(v);
                }
            }
        }
    }
    for (int i = 0; (unsigned) i < node->atoms.size(); ++i) {
        for (int j = 0; (unsigned) j < node->subgraphs.size(); ++j) {
            // searched for the current atom in the others
            if (node->subgraphs[j].contains(node->atoms[i])) {
                // if found, gets the ways to it
                way = node->subgraphs[j].get_paths_to(node->atoms[i]);
                for (auto& v : way) {
                    // adds the previous cut to every way returned by the
                    // get_paths_to function
                    v.insert(v.begin(), j);
                    // adds the way vector to the paths vector
                    paths.push_back(v);
                }
            }
        }
    }
}

std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    // 20p
    // creates the vector needed for the auxiliary function and
    // modifies it upon calling
    std::vector<std::vector<int>> paths;
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    // returns all the positions
    node->possible_deiterations_aux(paths);
    return paths;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // 10p
    AEGraph node_value = *this;
    AEGraph *node = &node_value;
    for (int i = 0; (unsigned) i < where.size() - 1; ++i) {
        node = &node->subgraphs[where[i]];
    }
    if (node->subgraphs.size() <= (unsigned) where[where.size() - 1]) {
        node->atoms.erase(node->atoms.begin() +
            (where[where.size() - 1] - node->subgraphs.size()));
    } else{
        node->subgraphs.erase(node->subgraphs.begin() +
            where[where.size() - 1]);
    }
    return node_value;
}

