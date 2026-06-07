class basriMSA:
    def __init__(self, sequences):
       
        self.sequences = sequences
        self.names = list(sequences.keys())
        self.num_seqs = len(sequences)
        self.distance_matrix = []
        self.guide_tree = None

    def calculate_distance_matrix(self):
        print("1. Mesafe Matrisi Hesaplanıyor (Needleman-Wunsch)...")
        matrix = [[0] * self.num_seqs for _ in range(self.num_seqs)]
        
        for i in range(self.num_seqs):
            for j in range(i + 1, self.num_seqs):
                seq1 = self.sequences[self.names[i]]
                seq2 = self.sequences[self.names[j]]
                
                distance = self.needleman_wunsch_distance(seq1, seq2)
                matrix[i][j] = distance
                matrix[j][i] = distance 
                
        self.distance_matrix = matrix
        return self.distance_matrix

    def needleman_wunsch_distance(self, seq1, seq2):
       
        match_cost = 0      
        mismatch_cost = 1   
        gap_penalty = 1     
        
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        
        for i in range(1, n + 1):
            dp[i][0] = dp[i-1][0] + gap_penalty
        for j in range(1, m + 1):
            dp[0][j] = dp[0][j-1] + gap_penalty
            
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                cost = match_cost if seq1[i-1] == seq2[j-1] else mismatch_cost
                
                match_mismatch = dp[i-1][j-1] + cost
                delete = dp[i-1][j] + gap_penalty
                insert = dp[i][j-1] + gap_penalty
                
                dp[i][j] = min(match_mismatch, delete, insert)
                
        return dp[n][m]

    def build_guide_tree(self):
        print("2. Kılavuz Ağaç Oluşturuluyor (UPGMA)...")
        clusters = {i: name for i, name in enumerate(self.names)}
        distance_dict = {}
        
        for i in range(self.num_seqs):
            for j in range(self.num_seqs):
                if i != j:
                    distance_dict[(i, j)] = self.distance_matrix[i][j]
                    
        active_clusters = list(clusters.keys())
        next_id = self.num_seqs
        
        while len(active_clusters) > 1:
            min_dist = float('inf')
            c1, c2 = -1, -1
            
            for i in range(len(active_clusters)):
                for j in range(i + 1, len(active_clusters)):
                    id1, id2 = active_clusters[i], active_clusters[j]
                    dist = distance_dict[(id1, id2)]
                    
                    if dist < min_dist:
                        min_dist, c1, c2 = dist, id1, id2
                        
            new_cluster_node = (clusters[c1], clusters[c2])
            clusters[next_id] = new_cluster_node
            
            for other_id in active_clusters:
                if other_id != c1 and other_id != c2:
                    new_dist = (distance_dict[(c1, other_id)] + distance_dict[(c2, other_id)]) / 2.0
                    distance_dict[(next_id, other_id)] = new_dist
                    distance_dict[(other_id, next_id)] = new_dist
                    
            active_clusters.remove(c1)
            active_clusters.remove(c2)
            active_clusters.append(next_id)
            next_id += 1
            
        self.guide_tree = clusters[active_clusters[0]]
        print(f"   Ağaç Hiyerarşisi: {self.guide_tree}")
        return self.guide_tree

    def progressive_alignment(self):
        print("3. Kademeli Hizalama (Progressive Alignment) Başlıyor...\n")
        return self._traverse_and_align(self.guide_tree)

    def _traverse_and_align(self, node):
        if isinstance(node, str):
            return {node: self.sequences[node]}
        
        left_node, right_node = node
        left_profile = self._traverse_and_align(left_node)
        right_profile = self._traverse_and_align(right_node)
        
        return self._align_profiles(left_profile, right_profile)

    def _get_consensus(self, profile):
        seqs = list(profile.values())
        if not seqs: return ""
        length = len(seqs[0])
        consensus = ""
        for i in range(length):
            col = [s[i] for s in seqs if s[i] != '-']
            if not col:
                consensus += '-'
            else:
                consensus += max(set(col), key=col.count)
        return consensus

    def _align_profiles(self, prof1, prof2):
        cons1 = self._get_consensus(prof1)
        cons2 = self._get_consensus(prof2)
        
        path = self._nw_traceback(cons1, cons2)
        aligned_profile = {}
        
        for name, seq in prof1.items():
            aligned_seq = ""
            idx = 0
            for op in path:
                if op in ['M', 'D']: 
                    aligned_seq += seq[idx] if idx < len(seq) else '-'
                    idx += 1
                else: 
                    aligned_seq += '-'
            aligned_profile[name] = aligned_seq
            
        for name, seq in prof2.items():
            aligned_seq = ""
            idx = 0
            for op in path:
                if op in ['M', 'I']: 
                    aligned_seq += seq[idx] if idx < len(seq) else '-'
                    idx += 1
                else: 
                    aligned_seq += '-'
            aligned_profile[name] = aligned_seq
            
        return aligned_profile

    def _nw_traceback(self, seq1, seq2):
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        tb = [[0] * (m + 1) for _ in range(n + 1)] 
        
        gap_pen = 1
        for i in range(1, n + 1):
            dp[i][0] = i * gap_pen
            tb[i][0] = 1 
        for j in range(1, m + 1):
            dp[0][j] = j * gap_pen
            tb[0][j] = 2 
            
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                cost = 0 if seq1[i-1] == seq2[j-1] else 1
                
                diag = dp[i-1][j-1] + cost
                up = dp[i-1][j] + gap_pen
                left = dp[i][j-1] + gap_pen
                
                min_val = min(diag, up, left)
                dp[i][j] = min_val
                
                if min_val == diag:
                    tb[i][j] = 0 
                elif min_val == up:
                    tb[i][j] = 1 
                else:
                    tb[i][j] = 2 
                    
        i, j = n, m
        path = []
        while i > 0 or j > 0:
            if tb[i][j] == 0:
                path.append('M') 
                i -= 1
                j -= 1
            elif tb[i][j] == 1:
                path.append('D') 
                i -= 1
            else:
                path.append('I') 
                j -= 1
                
        path.reverse()
        return path

    def run(self):
        self.calculate_distance_matrix()
        
        print("\nHesaplanan Mesafe Matrisi:")
        for row in self.distance_matrix:
            print(row)
        print("-" * 30)
        
        self.build_guide_tree()
        hizalanmis_diziler = self.progressive_alignment()
        
        print("Hizalama İşlemi Tamamlandı!")
        return hizalanmis_diziler