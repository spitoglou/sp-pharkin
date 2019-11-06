

def salt_factor(**kwargs):
    delivered_drug = kwargs.get('delivered_drug', False)
    dose_of_salt = kwargs.get('dose_of_salt', False)
    salt_factor = kwargs.get('salt_factor', False)
    if delivered_drug and salt_factor:
        return ('dose_of_salt', delivered_drug / salt_factor)
    if delivered_drug and dose_of_salt:
        return ('Salt Factor', delivered_drug / dose_of_salt)
    if dose_of_salt and salt_factor:
        return ('Delivered Drug', dose_of_salt * salt_factor)


if __name__ == '__main__':
    s = salt_factor(
        delivered_drug=400,
        dose_of_salt=500
    )
    print(s)
