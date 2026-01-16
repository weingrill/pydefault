from fpdf import FPDF


class PDF(FPDF):
    meine_seiten_zahl = 0

    def header(self):
        # Arial bold 15
        self.set_font('Arial', 'B', 15)
        # Move to the right
        self.cell(80)
        # Title
        self.cell(30, 10, 'Title', 1, 0, 'C')
        # Line break
        self.ln(20)

        # Page footer

    def footer(self):
        self.meine_seiten_zahl += 1
        # self.page_no()
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.meine_seiten_zahl) + '/{nb}', 0, 0, 'C')
        if self.meine_seiten_zahl >= 2:
            self.meine_seiten_zahl = 0
        # Instantiation of inherited class


pdf = PDF()
pdf.alias_nb_pages()
pdf.add_page()
pdf.set_font('Times', '', 12)
for i in range(1, 200):
    pdf.cell(0, 10, 'Printing line number ' + str(i), 0, 1)
    print(i, pdf.meine_seiten_zahl, pdf.page_no())
pdf.output('tuto2.pdf')
