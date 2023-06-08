import numpy as np

class GetVacancyFile:
    def __init__(self, vacancy_txt):
        self.vacancies = []
        self.total_vacancies = []
        with open(vacancy_txt, 'rU') as f:
            xy = []
            for j in f:
                m = j.rstrip().split('  ')

                if '-----------' in m:
                    break

            for j in f:
                m = j.rstrip().split('  ')

                if '' not in m:
                    self.vacancies.append([float(n) for n in m])
                else:
                    f.close()
                    break

            self.vacancies = np.array(self.vacancies)

            self.total_vacancies = np.zeros((len(self.vacancies[:,0]),2))
            self.total_vacancies[:,0] = self.vacancies[:,0]
            self.total_vacancies[:,1] = np.sum(self.vacancies[:,1:],axis=1)
