/*VI. Conditions of Use
 *
 * Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for Dent Earl and his institution to use the StatSigmMA-w software developed at the University of Washington for research purposes, on the following conditions:
 *
 * 1. Dent Earl and his institutional colleagues may modify the StatSigmMA-w software and distribute the resulting modified software for research purposes, provided (a) that it is distributed  together with these Conditions of Use, (b) that Martin Tompa receive a copy of the finalized modified software, and (c) that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with the authorship of the software.
 *  2. The StatSigmMA-w software will be used by you and/or your institution solely for noncommercial purposes, except with express permission from the authors.
 * 3. Any risk associated with using the StatSigmMA-w software at your institution is with you and your institution.
 * 4. StatSigmMA-w will be cited in any publication(s) reporting on data obtained from it as:
 *    Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome-Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
 *    Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567-572.
 */

#define d_MAX_LENGTH_NEWICK 4096
#define d_MAX_LENGTH_REF_SPECIES 32

int g_NUM_SPECIES = 0; // determined at runtime
class Tree;
Tree* root;
Tree** list_of_branches;
int branch_num;
int species_num;
char **species_name;

int name2branch(char* name) { 
    for (int i = 0; i < branch_num; i++)
        if ((strlen(species_name[i]) > 1) 
            && (strncmp(name, species_name[i], strlen(species_name[i])) == 0))
            return i;
    return -1;
}

class Tree {

 public:

    int left_num;
    int right_num;

    Tree* left_subtree;
    Tree* right_subtree;
    Tree* parent;

    char tree_name[1024];

    bool is_any_present;
    bool is_both_present;
    bool is_all_present;

    Tree(char* sp_name);
    Tree(Tree* l, Tree* r);
    ~Tree();
  
    void fill_present();
    void init_present();
};


void Tree::init_present() {
    is_any_present = is_both_present = is_all_present = false;
    if (left_subtree) {
        left_subtree->init_present();
        right_subtree->init_present();
    }
}


void Tree::fill_present() {
    if (left_subtree) {
        left_subtree->fill_present();
        right_subtree->fill_present();
        if (left_subtree->is_any_present && right_subtree->is_any_present)
            is_both_present = true;
        else
            is_both_present = false;

        if (left_subtree->is_any_present || right_subtree->is_any_present)
            is_any_present = true;
        else
            is_any_present = false;

        if (left_subtree->is_all_present && right_subtree->is_all_present)
            is_all_present = true;
        else
            is_all_present = false;
    }
}


Tree::Tree(char* sp_name) {
    left_subtree = NULL;
    right_subtree = NULL;
    parent = NULL;
    strcpy(tree_name, sp_name);
}


Tree::Tree(Tree* l, Tree* r) {
    left_subtree = l;
    l->parent = this;
    right_subtree = r;
    r->parent = this;

    parent = NULL;
    sprintf(tree_name,"(%s %s)",l->tree_name, r->tree_name);
}

Tree::~Tree() {
    if (left_subtree) delete left_subtree;
    if (right_subtree) delete right_subtree;
}

void checkPos(int& pos, unsigned limit)
{
    if ((int) limit < pos){
        cerr << "Error, pos > phylogeny length (" << pos 
             << " > " << limit << ")" << "." << endl
             << "This can be caused by a newick parsing error." << endl;
        exit(1);
    }
}

Tree* read_tree(char* newick_format, int& sp_num, int& pos) {
    while (isblank(newick_format[pos])){
        pos++;
        checkPos(pos, strlen(newick_format));
    }
    if (newick_format[pos] == '(') {
        pos++;
        checkPos(pos, strlen(newick_format));
        Tree* left = read_tree(newick_format, sp_num, pos);
        while (newick_format[pos] != ':'){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        pos++;
        checkPos(pos, strlen(newick_format));
        while (isblank(newick_format[pos])){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        int start_pos = pos;
        while (newick_format[pos] != ','){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        char branch_length[100];
        strncpy(branch_length, newick_format + start_pos, pos - start_pos);
        branch_length[pos - start_pos] = '\0';
        float left_branch = atof(branch_length);
        if (left_branch < 0.001) left_branch = 0.001;

        pos++;
        checkPos(pos, strlen(newick_format));
        Tree* right = read_tree(newick_format, sp_num, pos);
        while (newick_format[pos] != ':'){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        pos++;
        checkPos(pos, strlen(newick_format));
        while (isblank(newick_format[pos])){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        start_pos = pos;
        while (newick_format[pos] != ')'){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        strncpy(branch_length, newick_format + start_pos, pos - start_pos);
        branch_length[pos - start_pos] = '\0';
        float right_branch = atof(branch_length);
        if (right_branch < 0.001) right_branch = 0.001;
    
        pos++;
        checkPos(pos, strlen(newick_format));
        return (new Tree(left, right));
    }
    else {
        int start_pos = pos;
        while (newick_format[pos] != ':'){
            pos++;
            checkPos(pos, strlen(newick_format));
        }
        char sp_name[10000];
        strncpy(sp_name, newick_format + start_pos, pos - start_pos);
        sp_name[pos-start_pos] = '\0';    
        Tree* node = new Tree(sp_name);
        sp_num++;
        return node;
    }
}

void assign_branch_num(Tree* node) {
    Tree* ll = node->left_subtree;
    Tree* rr = node->right_subtree;
    
    if (node == root) {
        node->left_num = node->right_num = branch_num;
        list_of_branches[branch_num] = ll;
        // cout << branch_num << " " << ll->tree_name << endl;
        if (ll->left_subtree == NULL)
            strcpy(species_name[branch_num], ll->tree_name);
        if (rr->left_subtree == NULL)
            strcpy(species_name[branch_num], rr->tree_name);
        branch_num++;
    }
    else {
        node->left_num = branch_num;
        list_of_branches[branch_num] = ll;
        // cout << branch_num << " " << ll->tree_name << endl;
        if (ll->left_subtree == NULL)
            strcpy(species_name[branch_num], ll->tree_name);
        branch_num++;

        node->right_num = branch_num;
        list_of_branches[branch_num] = rr;
        // cout << branch_num << " " << rr->tree_name << endl;
        if (rr->left_subtree == NULL)
            strcpy(species_name[branch_num], rr->tree_name);
        branch_num++;
    }
    if (ll->left_subtree) assign_branch_num(ll);
    if (rr->left_subtree) assign_branch_num(rr);
}

unsigned countSpecies(Tree *node)
{
    // walks a tree and counts the number of leafs
    Tree* l = node->left_subtree;
    Tree* r = node->right_subtree;
    if ((l == NULL) && (r == NULL))
        return 1;
    return (countSpecies(l) + countSpecies(r));
}

int countNodes(Tree *node)
{
    // walks a tree and counts the number of nodes
    // based on the walking method present in assign_branch_num
    if (node == NULL)
        return 0;
    Tree* l = node->left_subtree;
    Tree* r = node->right_subtree;
    static int i = 0;
    int leftNum = 0, rightNum = 0, c = 0;
    if (node == root){
        if (g_options.printAll){
            cout << setw(2) << i << "\t" << l->tree_name << endl;
            ++i;
        }
        ++c;
    }else{
        if (g_options.printAll){
            cout << setw(2) << i << "\t" << l->tree_name << endl;
            ++i;
        }
        ++c;
        if (g_options.printAll){
            cout << setw(2) << i << "\t" << r->tree_name << endl;
            ++i;
        }
        ++c;
    }
    if (l->left_subtree) leftNum = countNodes(l);
    if (r->left_subtree) rightNum = countNodes(r);
    return (c + leftNum + rightNum);
}

void read_newick_tree(void) {
    int pos = 0;
    species_num = 0;
    root = read_tree(g_options.PHYLOGENY, species_num, pos);
    g_NUM_SPECIES = countNodes(root);
    list_of_branches = (Tree **) malloc(g_NUM_SPECIES * sizeof(Tree *));
    species_name = (char**) malloc(g_NUM_SPECIES * sizeof(char *));
    for (int i = 0; i < g_NUM_SPECIES; i++){
        species_name[i] = (char *) malloc(32);
        strcpy(species_name[i], "\0");
    }
    branch_num = 0;
    assign_branch_num(root);
}

bool is_OnPathToHuman(int leaf, int edge) {

    Tree* start = list_of_branches[leaf];
    Tree* ref = list_of_branches[name2branch(g_options.REF_SPECIES)];
    Tree* e = list_of_branches[edge];

    // Storing path from human to root

    Tree* path1[100];
    int len1 = 0;

    Tree* temp = ref;
    path1[len1++] = temp;
    while (temp->parent != NULL) {
        temp = temp->parent;
        path1[len1++] = temp;
    }


    // found denotes the common parent of that leaf and human
    // locate found from human's ancestors
    temp = start;
    int found=-1;
    for (int i = 0;i < len1;i++)
        if (temp == path1[i]) found = i;
  
    // locate found from leaf's ancestors,
    // also check whether edge is an ancestor of leaf
    do {
        // edge is an ancestor of leaf --> true
        if (temp == e) return true;
        temp = temp->parent;
        for (int i = 0;i < len1;i++)
            if (temp == path1[i]) found = i;
    } while (found==-1);

    // searching from human to found, check whether edge is on the path
    // ignoring the edge incident on the common parent (i.e. i = found is not considered)
    for (int i = 0;i < found;i++)
        if (path1[i] == e) return true;

    return false;
}

void insert_inner_branches(bool *cur_sp) {
    // Initialize
    root->init_present();
    // Set according to cur_sp
    for (int i = 0; i < branch_num; i++)
        list_of_branches[i]->is_both_present = list_of_branches[i]->is_any_present = cur_sp[i];

    Tree* ref = list_of_branches[name2branch(g_options.REF_SPECIES)];

    // Going up from human, and filling each sibling subtree
    Tree* temp = ref;
    Tree* prev;
    while (temp != root) {
        prev = temp->parent;
        if (prev->left_subtree == temp)
            prev->right_subtree->fill_present();
        else
            prev->left_subtree->fill_present();
        temp = prev;
    }

    // Now going up from reference species (i.e. human), and filling each on that path
    Tree* path[100];
    int len = 0;

    temp = ref;
    while (temp != NULL) {
        path[len++] = temp;
        temp = temp->parent;
    }

    // When judging if a internal branch is present, 
    // the refernce species is not counted.
    // For each node on the path from reference species to the root,
    // update is_both_present and is_any_present.
    // (Think the reference species as the root.)
    bool is_parent_present = false; // arbitrary initialization
    while (len > 1) {
        Tree* p = path[len - 1];

        // when p->left_subtree is on the path
        if (p->left_subtree == path[len - 2]) {
            // starting at root, where there is no sibling
            if (p == root) {
                // left subtree's properties is assigned by its sibling's properties.
                p->left_subtree->is_both_present = p->right_subtree->is_both_present;
                is_parent_present = p->left_subtree->is_any_present = p->right_subtree->is_any_present; 
            }
            // There is sibling
            else {
                // left subtree's properties is determined by those of its sibling and its ancestors' siblings.
                p->left_subtree->is_both_present = (is_parent_present && p->right_subtree->is_any_present);
                is_parent_present = p->left_subtree->is_any_present = (is_parent_present || p->right_subtree->is_any_present);
            }
        }

        // when p->right_subtree is on the path,
        // similar to p->left_subtree is on the path
        else {
            if (p == root) {
                p->right_subtree->is_both_present = p->left_subtree->is_both_present;
                is_parent_present = p->right_subtree->is_any_present =  p->left_subtree->is_any_present;
            }
            else {
                p->right_subtree->is_both_present = (is_parent_present && p->left_subtree->is_any_present);
                is_parent_present = p->right_subtree->is_any_present = (is_parent_present || p->left_subtree->is_any_present);
            }
        }
    
        len--;
    }

    // Only update from false -> true (Not from true -> false)
    // So it won't affect the branch leading to leaves (i.e. human).
    for (int s = 0; s < branch_num; s++)
        if (list_of_branches[s]->is_both_present)
            cur_sp[s] = true; 
}

void check_gapped_branches(bool *cur_sp) {
    // Initialize
    root->init_present();
    // Set according to cur_sp
    // is_all_present indicates if all leaves in the tree are present(i.e. gap).
    for (int s = 0; s < branch_num; s++)
        list_of_branches[s]->is_all_present = cur_sp[s];
    root->fill_present();
  
    // check if one subtree separated by current branch has all gaps at its leaves.
    for (int s = 0; s < branch_num; s++) {
        Tree* currTree = list_of_branches[s];
        bool allGaps = true;
        // If current tree (current branch leads to) is all gaps, then allGaps is true.
        // Otherwise, 
        if (!(currTree->is_all_present)) {
            Tree* temp = currTree;
            // check current tree's sibling and the sibling of each ancestor of current tree.
            while (temp != root) {
                Tree* prt = temp->parent;	
                if (temp == prt->left_subtree) {
                    if (!(prt->right_subtree->is_all_present)) {
                        allGaps = false;
                        break;
                    }
                }
                else {
                    if (!(prt->left_subtree->is_all_present)) {
                        allGaps = false;
                        break;
                    }
                }
                temp = prt;
            }
        }

        cur_sp[s] = allGaps;
    }
}
