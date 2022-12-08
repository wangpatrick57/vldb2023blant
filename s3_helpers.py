#!/bin/python3
class S3Manager:
    def __init__(self):
        self.s3_scores = dict()

    @staticmethod
    def get_species_in_order(speciesA, speciesB):
        return min(speciesA, speciesB), max(speciesA, speciesB)

    def add_s3_score(self, speciesA, speciesB, score):
        species1, species2 = S3Manager.get_species_in_order(speciesA, speciesB)

        if species1 not in self.s3_scores:
            self.s3_scores[species1] = dict()

        self.s3_scores[species1][species2] = score

    def get_s3_score(self, speciesA, speciesB):
        species1, species2 = S3Manager.get_species_in_order(speciesA, speciesB)
        return self.s3_scores[species1][species2]

    def add_all_s3_scores(self):
        self.add_s3_score('fly', 'chicken', 0.140055)
        self.add_s3_score('fly', 'duck', 0.151095)
        self.add_s3_score('fly', 'human', 0.152751)
        self.add_s3_score('fly', 'rat', 0.168291)
        self.add_s3_score('fly', 'cow', 0.171873)
        self.add_s3_score('fly', 'dog', 0.171989)
        self.add_s3_score('fly', 'guinea_pig', 0.179495)
        self.add_s3_score('fly', 'horse', 0.182391)
        self.add_s3_score('fly', 'cat', 0.183032)
        self.add_s3_score('fly', 'mouse', 0.186346)
        self.add_s3_score('chicken', 'rat', 0.560922)
        self.add_s3_score('rat', 'human', 0.599473)
        self.add_s3_score('cow', 'rat', 0.615501)
        self.add_s3_score('dog', 'rat', 0.637388)
        self.add_s3_score('cat', 'rat', 0.657103)
        self.add_s3_score('chicken', 'dog', 0.65752)
        self.add_s3_score('chicken', 'cow', 0.661516)
        self.add_s3_score('horse', 'rat', 0.664873)
        self.add_s3_score('duck', 'rat', 0.666449)
        self.add_s3_score('chicken', 'mouse', 0.67029)
        self.add_s3_score('chicken', 'cat', 0.671385)
        self.add_s3_score('chicken', 'horse', 0.676491)
        self.add_s3_score('guinea_pig', 'rat', 0.678282)
        self.add_s3_score('chicken', 'human', 0.714134)
        self.add_s3_score('cow', 'mouse', 0.714493)
        self.add_s3_score('dog', 'cow', 0.737947)
        self.add_s3_score('cat', 'cow', 0.759275)
        self.add_s3_score('cat', 'mouse', 0.776355)
        self.add_s3_score('guinea_pig', 'cow', 0.780501)
        self.add_s3_score('cat', 'dog', 0.781419)
        self.add_s3_score('guinea_pig', 'dog', 0.782205)
        self.add_s3_score('horse', 'dog', 0.78252)
        self.add_s3_score('cow', 'human', 0.782717)
        self.add_s3_score('duck', 'dog', 0.788338)
        self.add_s3_score('duck', 'cow', 0.789299)
        self.add_s3_score('duck', 'mouse', 0.794074)
        self.add_s3_score('duck', 'cat', 0.804088)
        self.add_s3_score('guinea_pig', 'cat', 0.804442)
        self.add_s3_score('guinea_pig', 'mouse', 0.806799)
        self.add_s3_score('cat', 'horse', 0.809356)
        self.add_s3_score('guinea_pig', 'horse', 0.824679)
        self.add_s3_score('dog', 'human', 0.831416)
        self.add_s3_score('duck', 'horse', 0.848655)
        self.add_s3_score('guinea_pig', 'human', 0.848896)
